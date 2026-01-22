"""
The purpose of this file is to populate for each field of a VCF line (and perform any modifications/calculations to achieve this)
"""

from Bio import SeqIO
from convert_gvf_to_vcf.assistingconverter import convert_gvf_attributes_to_vcf_values
from convert_gvf_to_vcf.logger import logger
from dataclasses import dataclass
from typing import Optional,Union


def extract_reference_allele(fasta_file, chromosome_name, position, end):
    """ Extracts the reference allele from the assembly.
    :param fasta_file: FASTA file of the assembly
    :param chromosome_name: name of the sequence
    :param position: position
    :param end: end position
    :return: reference_allele: base found at this chromosome_name at this position within this fasta_file
    """
    # using .index for memory efficiency:  https://biopython.org/docs/1.76/api/Bio.SeqIO.html#input-multiple-records
    records_dictionary = SeqIO.index(fasta_file, "fasta")
    zero_indexed_position = position - 1  # minus one because zero indexed
    zero_indexed_end = end - 1
    reference_allele = ""
    for position in range(zero_indexed_position, zero_indexed_end):
        reference_allele = reference_allele + records_dictionary[chromosome_name].seq[position]
    records_dictionary.close()
    return reference_allele

@dataclass
class VariantRange:
    """
    This is responsible for representing the co-ordinate ranges of a structural variant.
    Attributes:
        pos: start
        start_range_lower_bound: can be int "." or None
        start_range_upper_bound: can be int "." or None
        end: end
        end_range_lower_bound: can be int "." or None
        end_range_upper_bound: can be int "." or None
    """
    pos: int
    start_range_lower_bound: Optional[Union[int, str]] #value can be int, "." or None
    start_range_upper_bound: Optional[Union[int, str]] #value can be int, "." or None
    end: int
    end_range_lower_bound: Optional[Union[int, str]] #value can be int, "." or None
    end_range_upper_bound: Optional[Union[int, str]] #value can be int, "." or None

class VcfLineBuilder:
    """
    This class is responsible for creating VcfLine objects that contain VCF datalines.
    """
    def __init__(self, field_lines_dictionary, all_possible_lines_dictionary, reference_lookup, ordered_list_of_samples):
        self.field_lines_dictionary = field_lines_dictionary
        self.all_possible_lines_dictionary = all_possible_lines_dictionary
        self.reference_lookup = reference_lookup
        self.ordered_list_of_samples = ordered_list_of_samples

    def build_vcf_line(self, gvf_feature_line_object):
        # Attributes which store important key-values dicts
        (vcf_value_from_gvf_attribute,  # used to populate the VCF fields. This is a dict of non-converted GVF attribute keys and their values.
         vcf_values_for_info,  # a dict that stores INFO key-values to form VCF line. This includes converted GVF attribute keys (+ other SV INFO).
         vcf_values_for_format  # a dict of FORMAT key-values for each sample to form VCF line
         ) = convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes, self.reference_lookup.mapping_attribute_dict, self.field_lines_dictionary, self.all_possible_lines_dictionary)

        # Attributes which might form useful parts of INFO field in VCF lines (useful information from GVF)
        source = gvf_feature_line_object.source
        so_type = gvf_feature_line_object.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?
        end = int(gvf_feature_line_object.end)
        phase = gvf_feature_line_object.phase # this is always a placeholder '.'

        # Attributes which are required to generate a VCF DATALINE
        # MANDATORY VCF FIELD 1
        chrom = gvf_feature_line_object.seqid
        # MANDATORY VCF FIELD 2
        pos = int(gvf_feature_line_object.start)
        # MANDATORY VCF FIELD 3
        id = vcf_value_from_gvf_attribute["ID"] # attributes: ID
        # note ref and alt are calculated below (fields 4 and 5)
        # MANDATORY VCF FIELD 6
        qual = gvf_feature_line_object.score  # see EVA-3879: this is always '.'
        # MANDATORY VCF FIELD 7
        filter = "." # this is always a placeholder '.'
        # forms MANDATORY VCF FIELD 8
        info_dict = {}  # dict that stores all INFO key-values (including INFO from merged lines and SV INFO).
        # calculated last
        length = end - pos # required for INFO fields- SVLEN and END
        # MANDATORY VCF FIELD 4
        ref = self.get_ref(vcf_value_from_gvf_attribute, chrom, pos, end)
        # MANDATORY VCF FIELD 5
        pos, ref, alt, info_dict = self.get_alt(vcf_value_from_gvf_attribute, chrom, pos, end, length, ref, so_type)
        if info_dict:
            vcf_values_for_info.update(info_dict)
        return VcfLine(chrom=chrom, pos=pos, id=id, ref=ref, alt=alt, qual=qual, filter=filter,
                       info_dict=vcf_values_for_info, vcf_values_for_format=vcf_values_for_format,
                       order_sample_names=self.ordered_list_of_samples)

    def build_vcf_header(self):
        """Stores VCF header lines and makes them unique.
        :return" self.field_lines_dictionary: VCF header lines ALT,INFO,FILTER,FORMAT
        """
        # make unique
        self.field_lines_dictionary["ALT"] = list(dict.fromkeys(self.field_lines_dictionary["ALT"]))
        self.field_lines_dictionary["INFO"] = list(dict.fromkeys(self.field_lines_dictionary["INFO"]))
        self.field_lines_dictionary["FILTER"] = list(dict.fromkeys(self.field_lines_dictionary["FILTER"]))
        self.field_lines_dictionary["FORMAT"] = list(dict.fromkeys(self.field_lines_dictionary["FORMAT"]))
        return self.field_lines_dictionary

    # Functions which are responsible for token generation/population for the VCF line
    def add_padded_base(self, chrom, pos, end, ref, alt, placed_before : bool):
        """ Adds a padded base to the REF and ALT allele of a VCF line.
        :param ref: reference allele
        :param alt: alt allele
        :param placed_before: padded base is placed before ref or alt True or False
        :return: (padded_base, pos, ref, alt)
        """
        if placed_before:
            pos = pos - 1
            padded_base = extract_reference_allele(self.reference_lookup.assembly_file, chrom, pos, pos + 1)
            ref = padded_base + ref
        else:
            end = end + 1
            padded_base = extract_reference_allele(self.reference_lookup.assembly_file, chrom, end-1, end)
            ref = ref + padded_base
        return padded_base, pos, ref, alt

    def convert_iupac_ambiguity_code(self, ref_to_convert):
        """ If the REF allele of a VCF line contains an IUPAC ambiguity code, converts it.
        :param ref_to_convert: reference allele to be converted
        :return: self.ref
        """
        converted_ref = ""
        for base in ref_to_convert:
            if base in self.reference_lookup.iupac_ambiguity_dictionary:
                iupac_value = min(self.reference_lookup.iupac_ambiguity_dictionary[base])
                converted_base = iupac_value
            else:
                converted_base = base
            converted_ref = converted_ref + converted_base
        return converted_ref

    def check_ref(self, ref_allele_to_be_checked):
        """ Checks whether a reference allele meets the requirements of the VCF specification.
        :param ref_allele_to_be_checked: reference allele to check
        :return: checked_reference_allele: reference allele that meets the requirements of the VCF specification"""
        if isinstance(ref_allele_to_be_checked, str):
            if not all(bases in ref_allele_to_be_checked for bases in ["A", "C", "G", "T", "N"]):
                checked_reference_allele = self.convert_iupac_ambiguity_code(ref_allele_to_be_checked)
            else:
                checked_reference_allele = ref_allele_to_be_checked
        else:
            logger.warning(f"Ref allele must be a string: {ref_allele_to_be_checked}. Setting to missing value.")
            checked_reference_allele = "."
        return checked_reference_allele


    def get_ref(self, vcf_value_from_gvf_attribute, chrom, pos, end):
        """ Gets the reference allele from attributes column or if not found, returns "."
        :return: reference allele
        """
        assembly_file = self.reference_lookup.assembly_file
        if "Reference_seq" in vcf_value_from_gvf_attribute.keys():
            reference_allele = vcf_value_from_gvf_attribute["Reference_seq"]
        else:
            if assembly_file:
                reference_allele = extract_reference_allele(assembly_file, chrom, pos, end)
            else:
                logger.warning("No reference provided. Placeholder inserted for Reference allele.")
                reference_allele = "."
        if reference_allele != ".":
            reference_allele = self.check_ref(reference_allele)
        return reference_allele

    def create_coordinate_range(self, vcf_value_from_gvf_attribute):
        """ Create the start and end range using the dictionary of GVF attributes and the pos and end
        :return: start_range_lower_bound, start_range_upper_bound, end_range_lower_bound, end_range_upper_bound
        """
        if "Start_range" in vcf_value_from_gvf_attribute:
            # Set value to the start range found in the GVF attribute
            start_range_lower_bound = vcf_value_from_gvf_attribute["Start_range"][0]
            start_range_upper_bound = vcf_value_from_gvf_attribute["Start_range"][1]
        else:
            # Setting values for a precise variant. Start Ranges only apply to imprecise variants.
            start_range_lower_bound = None
            start_range_upper_bound = None
        if "End_range" in vcf_value_from_gvf_attribute:
            # Set value to the end range found in the GVF attribute
            end_range_lower_bound = vcf_value_from_gvf_attribute["End_range"][0]
            end_range_upper_bound = vcf_value_from_gvf_attribute["End_range"][1]
        else:
            # Set value for a precise variant. End Ranges only apply to imprecise variants.
            end_range_lower_bound = None
            end_range_upper_bound = None
        return start_range_lower_bound, start_range_upper_bound, end_range_lower_bound, end_range_upper_bound

    def generate_info_field_symbolic_allele(self,
                                            variant_range_coordinates,
                                            length, ref, symbolic_allele):
        """ Generate the dictionary INFO field values for the symbolic allele.
        :param variant_range_coordinates: object with attribues:pos, start_range_lower_bound, start_range_upper_bound, end, end_range_lower_bound, end_range_upper_bound
        :param length:  end - pos
        :param ref: reference allele
        :param symbolic_allele: symbolic allele for ALT
        :return:
            info_dict: dict of key-value pairs for INFO field (END, IMPRECISE, CIPOS, CIEND, SVLEN),
            is_imprecise: boolean
        """

        # Set up info dictionary
        info_dict = {}
        # Setting up the info keys with None as a default value
        for info_key in ["END", "IMPRECISE", "CIPOS", "CIEND", "SVLEN"]:
            info_dict[info_key] = None
        if length:
            info_dict["SVLEN"] = str(length)
        # Calculate precise or imprecise
        range_keys = ["start_range_lower_bound", "start_range_upper_bound", "end_range_lower_bound", "end_range_upper_bound"]
        # in the case of a PRECISE variant
        if any(range_key is None for range_key in range_keys):
            info_dict["END"], is_imprecise = self.generate_info_field_for_precise_variant(variant_range_coordinates.pos, ref)
        # IMPRECISE variants
        else:
            info_fields_for_imprecise_variants = self.generate_info_field_for_imprecise_variant(
                variant_range_coordinates,
                length, ref,
                symbolic_allele
            )
            info_dict.update({
                "CIEND": info_fields_for_imprecise_variants[0],
                "CIPOS": info_fields_for_imprecise_variants[1],
                "END": info_fields_for_imprecise_variants[2],
                "IMPRECISE": info_fields_for_imprecise_variants[3]})
            # Set the imprecise flag
            is_imprecise = info_fields_for_imprecise_variants[4]
        return info_dict, is_imprecise

    def generate_info_field_for_imprecise_variant(self,
                                                  variant_range_coordinates,
                                                  length, ref,
                                                  symbolic_allele):
        """ Generate the INFO field values for an imprecise variant.
        :param variant_range_coordinates: object with attribues:pos, start_range_lower_bound, start_range_upper_bound, end, end_range_lower_bound, end_range_upper_bound
        :param length:  end - pos
        :param ref: refernece allele
        :param symbolic_allele: symbolic allele for ALT
        :return:
            info_ciend_value(string i.e "0,100"),
            info_cipos_value(string i.e "0,1"),
            info_end_value (string i.e. "13"),
            info_imprecise_value(string=="IMPRECISE"),
            is_imprecise (boolean == True)
        """
        # Imprecise variant
        is_imprecise = True
        info_imprecise_value = "IMPRECISE"
        ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound = self.calculate_CIPOS_and_CIEND(
            variant_range_coordinates)
            # end, end_range_lower_bound, end_range_upper_bound, pos, start_range_lower_bound, start_range_upper_bound)

        # form the CIPOS value
        info_cipos_value =  f"{str(cipos_lower_bound)},{str(cipos_upper_bound)}" if all(cipos_bound is not None for cipos_bound in (cipos_lower_bound, cipos_upper_bound)) else None
        # form the CIEND value
        info_ciend_value = f"{str(ciend_lower_bound)},{str(ciend_upper_bound)}" if all(ciend_bound is not None for ciend_bound in (cipos_lower_bound, cipos_upper_bound)) else None

        # Determine the END value based on the symbolic allele
        if symbolic_allele == "<INS>":
            info_end_value = str(variant_range_coordinates.pos + len(ref) - 1)
        elif any(sym in symbolic_allele for sym in {"DEL", "DUP", "INV", "CNV"}):
            info_end_value = str(variant_range_coordinates.pos + length)
        elif symbolic_allele == "<*>":
            info_end_value = str(variant_range_coordinates.pos + len(ref))
        else:
            logger.warning(f"Cannot identify symbolic allele: {symbolic_allele}")
        return info_ciend_value, info_cipos_value, info_end_value, info_imprecise_value, is_imprecise

    def generate_info_field_for_precise_variant(self, pos, ref):
        """ Creates values for the VCF file INFO field if the variant is precise
        :param pos: gvf_feature_line_object.start
        :param ref: reference allele
        :return: info_end_value(string), is_imprecise(boolean)
        """
        is_imprecise = False
        info_end_value = str(pos + len(ref) - 1)
        return info_end_value, is_imprecise

    def calculate_CIPOS_and_CIEND(self,
                                  variant_range_cordinates):
        """Converts start and end range co-ordinates into CIPOS and CIEND (for VCF info field).
        :param variant_range_cordinates:  object with attributes:pos, start_range_lower_bound, start_range_upper_bound, end, end_range_lower_bound, end_range_upper_bound
        :return: ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound
        """
        # Converting start_range co-ordinates
        cipos_lower_bound = self.convert_to_ci_bound(variant_range_cordinates.pos, variant_range_cordinates.start_range_lower_bound)
        cipos_upper_bound = self.convert_to_ci_bound(variant_range_cordinates.pos, variant_range_cordinates.start_range_upper_bound)
        # Converting End_range co-ordinates
        ciend_lower_bound = self.convert_to_ci_bound(variant_range_cordinates.end, variant_range_cordinates.end_range_lower_bound)
        ciend_upper_bound = self.convert_to_ci_bound(variant_range_cordinates.end, variant_range_cordinates.end_range_upper_bound)
        return ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound

    def convert_to_ci_bound(self, start_or_end_position, original_coordinate_bound):
        """ Logic for converting start or end range co-ordinates to CIPOS or CIEND
        :param: start or end positon: (pos or end)
        :param: original_coordinate bound: start/end range lower or upper bound
        :param: cipos/ciend upper or lower bound
        """
        # if the original bound is unknown for the imprecise variant, pass this to CIPOS or CIEND(i.e. CIPOS=.,0)
        if original_coordinate_bound == ".":
            ci_coordinate_bound = "."
        # if the original bound is precise variant, set to CIPOS or CIEND to None
        elif original_coordinate_bound is None:
            ci_coordinate_bound = None
        # if known imprecise coordinate, calculate the CIPOS/CIEND
        else:
            ci_coordinate_bound = int(original_coordinate_bound) - start_or_end_position
        return ci_coordinate_bound

    def generate_symbolic_allele(self, vcf_value_from_gvf_attribute, pos, end, length, ref, so_type):
        """ Generates the symbolic allele and stores the corresponding metainformation lines.
        Also determines if variant is precise or imprecise.
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        # Get symbolic allele
        symbolic_allele_id = self.reference_lookup.symbolic_allele_dictionary[so_type][1]
        symbolic_allele = f'<{symbolic_allele_id}>'

        lines_standard_alt = self.field_lines_dictionary["ALT"]
        lines_standard_info = self.field_lines_dictionary["INFO"]
        all_possible_alt_lines = self.all_possible_lines_dictionary["ALT"]
        all_possible_info_lines = self.all_possible_lines_dictionary["INFO"]

        if symbolic_allele_id in all_possible_alt_lines:
            lines_standard_alt.append(all_possible_alt_lines[symbolic_allele_id])

        # Creating start/end co-ordinate ranges
        (start_range_lower_bound, start_range_upper_bound,
         end_range_lower_bound, end_range_upper_bound) = self.create_coordinate_range(vcf_value_from_gvf_attribute)
        variant_range_coordinates = VariantRange(pos=pos,
                                                 start_range_lower_bound=start_range_lower_bound,
                                                 start_range_upper_bound=start_range_upper_bound,
                                                 end=end,
                                                 end_range_lower_bound=end_range_lower_bound,
                                                 end_range_upper_bound=end_range_upper_bound)

        # Get the dictionary of INFO keys and their values
        info_dict, is_imprecise = self.generate_info_field_symbolic_allele(
            variant_range_coordinates, length, ref, symbolic_allele)

        # Set SVCLAIM depending on the symbolic allele
        info_svclaim_key = "SVCLAIM"
        if "DEL" in symbolic_allele or "DUP" in symbolic_allele:
            # TODO: IMPORTANT: this should be set to the missing place holder '.' but await clarification from the spec
            info_svclaim_value = "D"
            info_dict.update({info_svclaim_key: info_svclaim_value})
            lines_standard_info.append(all_possible_info_lines["SVCLAIM"])

        # for all variants (precise and imprecise) store INFO lines for the header
        lines_standard_info.append(all_possible_info_lines["END"])
        lines_standard_info.append(all_possible_info_lines["SVLEN"])

        # for imprecise variants only
        if is_imprecise:
            lines_standard_info.append(all_possible_info_lines["IMPRECISE"])
            lines_standard_info.append(all_possible_info_lines["CIPOS"])
            lines_standard_info.append(all_possible_info_lines["CIEND"])
        return symbolic_allele, info_dict, lines_standard_alt, lines_standard_info

    #TODO: awaiting clarification from the specification to keep or remove this function
    def has_svclaim_abundance_evidence(self, vcf_value_from_gvf_attribute, alt, info_dict):
        """ This functions determines if the GVF attributes contain evidence that the SVCLAIM should be set to Abundance (D)
        :params: vcf_value_from_gvf_attribute: dictionary of GVF attributes
        :params: alt: alternative allele (SVCLAIM mandatory for DEL/DUP)
        :params: info_dict: dictionary of INFO keys and their values
        return: boolean/none
        """
        abundance_variant_method_keywords = ["Array_CGH", "Read_depth", "Depth_of_coverage", "SNP_Microarray", "Oligo_Microarray"]
        abundance_variant_method_keywords_lowercase = {vm.casefold() for vm in abundance_variant_method_keywords}

        abundance_variant_region_description_keywords = ["Inferred", "Copy_number_variation", "CNV", "Gain", "Loss", "Dosage_sensitive", "Relative_copy_number"]
        abundance_variant_region_description_keywords_lowercase = {vrd.casefold() for vrd in abundance_variant_region_description_keywords}

        abundance_variant_call_description_keywords = ["Inferred", "Copy_number_variation", "CNV", "Gain", "Loss", "Dosage_sensitive", "Relative_copy_number"]
        abundance_variant_call_description_keywords_lowercase = {vcd.casefold() for vcd in abundance_variant_call_description_keywords}

        abundance_variant_seq_keywords = ["."]
        # Region
        # SO:0001019 = copy_number_variation
        abundance_variant_region_so_id_keywords = ["SO:0001019"]

        # Call
        # SO:0001742 = copy_number_gain
        # SO:0001743 = copy_number_loss
        # SO:1000173 = tandem_duplication
        abundance_variant_call_so_id_keywords =   ["SO:0001742", "SO:0001743", "SO:1000173"]
        abundance_imprecise_keywords = [True]

        if "DEL" in alt or "DUP" in alt:
            has_abundance_variant_method = vcf_value_from_gvf_attribute.get("Variant_Method", "").casefold() in abundance_variant_method_keywords_lowercase
            has_abundance_variant_region_description = any(i in vcf_value_from_gvf_attribute.get("variant_region_description", "").casefold() for i in abundance_variant_region_description_keywords_lowercase)
            has_abundance_call_description = any(i in vcf_value_from_gvf_attribute.get("variant_call_description", "").casefold() for i in abundance_variant_call_description_keywords_lowercase)
            has_abundance_variant_seq = vcf_value_from_gvf_attribute.get("Variant_seq") in abundance_variant_seq_keywords
            has_abundance_variant_region_so_id_keyword = vcf_value_from_gvf_attribute.get("variant_region_so_id") in abundance_variant_region_so_id_keywords
            has_abundance_variant_call_so_id_keyword = vcf_value_from_gvf_attribute.get("variant_call_so_id") in abundance_variant_call_so_id_keywords
            has_abundance_imprecise = info_dict.get("IMPRECISE") in abundance_imprecise_keywords
            if any([has_abundance_variant_method, has_abundance_variant_region_description,has_abundance_call_description, has_abundance_variant_seq, has_abundance_variant_region_so_id_keyword, has_abundance_variant_call_so_id_keyword, has_abundance_imprecise]):
                return True
            else:
                return False
        return None

    def get_alt(self, vcf_value_from_gvf_attribute, chrom, pos, end, length, ref, so_type):
        """ Gets the ALT allele for the VCF file
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        info_dict = {}
        if any(base in vcf_value_from_gvf_attribute["Variant_seq"] for base in ["A", "C", "G", "T", "N"]):
            alt = vcf_value_from_gvf_attribute["Variant_seq"]
        elif vcf_value_from_gvf_attribute["Variant_seq"] == '.':
            symbolic_allele, info_dict, lines_standard_alt, lines_standard_info = self.generate_symbolic_allele(vcf_value_from_gvf_attribute, pos, end, length, ref, so_type)

            if symbolic_allele is None:
                alt = "."
            elif (vcf_value_from_gvf_attribute["Variant_seq"] == "." or vcf_value_from_gvf_attribute["Variant_seq"] == "-") and symbolic_allele is not None:
                alt = symbolic_allele
                # add padded bases
                if pos == 1:
                    #print("pos, ref, alt",self.pos,self.ref, alt)
                    padded_base, pos, ref, alt = self.add_padded_base(chrom, pos, end, ref, alt, False)
                    ref = self.check_ref(ref)
                else:
                    #print("pos, ref, alt", self.pos,self.ref, alt)
                    padded_base, pos, ref, alt = self.add_padded_base(chrom, pos, end, ref, alt, True)
                    ref = self.check_ref(ref)
            else:
                alt = "."
                logger.warning(f"Cannot identify symbolic allele: {symbolic_allele}. Variant type is not supported.")
        else:
            alt = "."
            logger.warning("Could not determine the alternative allele.")
        #TODO: awaiting clarification to keep or remove this for SVCLAIM
        #is_abundance = self.has_svclaim_abundance_evidence(vcf_value_from_gvf_attribute, alt, info_dict)
        return pos, ref, alt, info_dict


class VcfLine:

    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info_dict=None, vcf_values_for_format=None,
                 order_info_list=None, order_format_keys=None, order_sample_names=None):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        if info_dict:
            self.info_dict = info_dict
        else:
            self.info_dict = {}
        if vcf_values_for_format:
            self.vcf_values_for_format = vcf_values_for_format
        else:
            self.vcf_values_for_format = {}

        # TODO: Use this to enforce the order of the info key values
        self.order_info_list = order_info_list if order_info_list else list(self.info_dict.keys())
        self.order_sample_names = order_sample_names if order_sample_names  else list(self.vcf_values_for_format.keys())
        self.order_format_keys = order_format_keys


    def __str__(self):
        """ Creates and formats the VCF line.
        :return: string_to_return - the VCF line as a string
        """
        fields = [self.chrom,
                  str(self.pos),
                  self.id,
                  self.ref,
                  self.alt,
                  self.qual,
                  self.filter,
                  self.format_info_string()]
        if self.format_keys != ["."]:
            fields.append(":".join(self.format_keys))
            fields.append(self.combine_format_values_by_sample_as_str())
        string_to_return = '\t'.join(fields)
        return string_to_return

    def __eq__(self, other_vcf_line):
        """ Compares equality of PARTS of the VcfLine objects.
        :param: other_vcf_line: another object to compare equality with
        """
        if isinstance(other_vcf_line, VcfLine):
            return (self.chrom == other_vcf_line.chrom) and (self.pos == other_vcf_line.pos) and (self.ref == other_vcf_line.ref)
        return False

    @property
    def format_keys(self):
        set_of_format_key = set([format_key
                                 for format_value in self.vcf_values_for_format.values()
                                 for format_key in format_value.keys()
                                 ])
        if len(set_of_format_key) == 0:
            set_of_format_key = set('.')
        return self._order_format_keys(set_of_format_key)  # a list of ordered format keys

    def merge_and_add(self, previous_element, current_element, delimiter):
        """ Merges fields of a VCF line. If field is the same, use current element. If different, merge with delimiter.
        :param: previous_element
        :param: current_element
        :param: delimiter
        :return: merged element
        """
        if previous_element == current_element:
            merged_element = current_element
        else:
            merged_element = delimiter.join((previous_element, current_element))
        return merged_element

    # functions responsible for FORMAT are below
    def _order_format_keys(self, set_of_format_keys):
        """Stores the FORMAT keys of the VCF line in the correct order by anchoring GT as the first key.
        :param: set_of_format_keys: format keys in a set
        :return: anchored_list_of_keys: list of ordered keys
        """
        if self.order_format_keys:
            return self.order_format_keys
        anchored_list_of_format_keys = []
        if 'GT' in set_of_format_keys:
            anchored_list_of_format_keys.append("GT")
            set_of_format_keys.remove('GT')
        anchored_list_of_format_keys.extend(sorted(set_of_format_keys))
        return anchored_list_of_format_keys

    def combine_format_values_by_sample_as_str(self):
        """ Creates list of format values for each sample for the vcf data line.
        :return: formatted string containing each format column e.g. '.:3\t.:.\t.:.\t0:1:3' (in the VCF file, this would be the tab-separated values under the sample name)
        """
        format_keys = self.format_keys
        list_of_sample_names = self.order_sample_names
        list_of_format_values_per_sample = []
        # Generate string. For present samples, get its format value. For missing samples, populate with a missing value.
        for sample in list_of_sample_names:
            if sample in self.vcf_values_for_format:
                format_value_list = []
                for key in format_keys:
                    format_value_list.append(
                        str(self.vcf_values_for_format.get(sample,'.').get(key,'.')) # adds missing values if not found
                    )
                list_of_format_values_per_sample.append(":".join(format_value_list))
            else:
                list_of_format_values_per_sample.append(':'.join(['.' for key in format_keys] or ['.']))
        return '\t'.join(list_of_format_values_per_sample)

    # functions responsible for INFO are below

    def fill_merge_dicts(self, merged_info_dict, key, previous_line_info_value, current_line_info_value):
        """ Logic for merging info dicts
        :param: merged_info_dict: merged dictionary
        :param: key: key for merged info dict
        :param: previous_line_info_value
        :param: current_line_info_value
        :return: merged info dict
        """
        if previous_line_info_value is None and current_line_info_value is None:
            pass
        elif previous_line_info_value == current_line_info_value:
            merged_info_dict[key] = previous_line_info_value
        else:
            if previous_line_info_value is None:
                merged_info_dict[key] = current_line_info_value
            elif current_line_info_value is None:
                merged_info_dict[key] = previous_line_info_value
            else:
                merged_info_dict[key] = f"{previous_line_info_value},{current_line_info_value}"
        return merged_info_dict

    def merge_info_dicts(self, other_vcf_line):
        """ Merges and stores the INFO dictionaries for the INFO field of a VCF line.
        :param: other_vcf_line
        :return: None
        """
        # Create data structure to merge the INFO dict of this VCF line and the other_vcf_line
        merged_info_dict = {}
        # Step 1: Merge from info_dict
        # Aim is to store SV INFO
        # info_dict = dict that stores all INFO key-values (including INFO from merged lines and SV INFO).
        for info_dict_key in self.info_dict.keys() | other_vcf_line.info_dict.keys():
            this_info_dict_value = self.info_dict.get(info_dict_key)
            other_info_dict_value = other_vcf_line.info_dict.get(info_dict_key)
            merged_info_dict = self.fill_merge_dicts(merged_info_dict,info_dict_key, this_info_dict_value,other_info_dict_value)

        # Remove the ID
        key_to_remove = "ID"
        if key_to_remove in merged_info_dict:
            del merged_info_dict[key_to_remove]

        # Store merged info dict for this VCF line and the other VCF line.
        self.info_dict = merged_info_dict
        other_vcf_line.info_dict = merged_info_dict

    def merge_vcf_values_for_format(self, other_vcf_line):
        """ Merge FORMAT vcf value dictionary for two vcf lines.
        :param: other_vcf_line
        return: None
        """
        for sample in other_vcf_line.vcf_values_for_format:
            if sample not in self.vcf_values_for_format:
                self.vcf_values_for_format[sample] = other_vcf_line.vcf_values_for_format.get(sample, {})
            else:
                for format_key in other_vcf_line.vcf_values_for_format.get(sample):
                    if format_key not in self.vcf_values_for_format[sample]:
                        self.vcf_values_for_format[sample][format_key] = other_vcf_line.vcf_values_for_format[sample][format_key]

    def format_info_string(self):
        """ Creates a formatted INFO string using the INFO dictionary. Anchors ID to start of the string.
        :return: info_string: formatted INFO string for use in VCF line
        """
        # Ensure ID is the first key
        anchored_key = "ID"
        if anchored_key in self.info_dict:
            self.info_dict = {anchored_key:self.info_dict.pop(anchored_key), **self.info_dict}
        # Remove None values
        if "IMPRECISE" in self.info_dict and self.info_dict.get('IMPRECISE') is None:
            del self.info_dict["IMPRECISE"]
        if "CIPOS" in self.info_dict and self.info_dict.get('CIPOS') is None:
            del self.info_dict["CIPOS"]
        if "CIEND" in self.info_dict and self.info_dict.get('CIEND') is None:
            del self.info_dict["CIEND"]
        # Format the string
        info_string = ";".join(f"{key}={value}" if key != "IMPRECISE" else f"{value}" for key,value in self.info_dict.items())
        return info_string

    # MERGE OR KEEP below
    def merge(self, other_vcf_line, list_of_sample_names):
        """ Merging of the fields of a VCF line (ID, ALT, FILTER, INFO, FORMAT, FORMATvalues).
        :param: other_vcf_line : other VCF line to merge with
        :param: list_of_sample_names: list of sample names to help with creating format values by sample
        """
        # Merging ID, ALT and FILTER first
        merged_id = self.merge_and_add(self.id, other_vcf_line.id, ";")
        sorted_merged_id_set = sorted(set(map(int, merged_id.split(";"))))
        sorted_merged_id = ";".join(map(str, sorted_merged_id_set))
        merged_alt = self.merge_and_add(self.alt, other_vcf_line.alt, ",")
        merged_filter = self.merge_and_add(self.filter, other_vcf_line.filter, ";")

        self.id = other_vcf_line.id = sorted_merged_id
        self.alt = other_vcf_line.alt = merged_alt
        self.filter = other_vcf_line.filter = merged_filter
        # Merging INFO using info_dict
        self.merge_info_dicts(other_vcf_line)
        # Merging FORMAT values - these go under the Sample
        self.merge_vcf_values_for_format(other_vcf_line)
        return self


    def keep(self, list_of_sample_names):
        self.order_sample_names = list_of_sample_names
        return self
