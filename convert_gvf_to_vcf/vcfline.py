# the purpose of this file is to populate for each field of a VCF line (and perform any necessary calculations to achieve this)
from convert_gvf_to_vcf.assistingconverter import convert_gvf_attributes_to_vcf_values
from Bio import SeqIO

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

class VcfLine:
    def __init__(self, gvf_feature_line_object,
                 info_attribute_dict,
                 symbolic_allele_dictionary,
                 assembly_file,
                 field_lines_dictionary,
                 all_possible_lines_dictionary):
        # assisting_converter = Assistingconverter()
        # self.vcf_value, self.info_string = assisting_converter.convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes, info_attribute_dict, field_lines_dictionary, all_possible_lines_dictionary)
        self.vcf_value, self.info_string = convert_gvf_attributes_to_vcf_values(gvf_feature_line_object.attributes, info_attribute_dict, field_lines_dictionary, all_possible_lines_dictionary)
        # ATTRIBUTES
        self.assembly = assembly_file
        self.symbolic_allele_dictionary = symbolic_allele_dictionary
        self.iupac_ambiguity_dictionary = self.build_iupac_ambiguity_code()
        # GVF
        self.source = gvf_feature_line_object.source
        self.so_type = gvf_feature_line_object.feature_type #currently column 3 of gvf, but could be an attribute so perhapsVCF: INFO or FORMAT?
        self.end = int(gvf_feature_line_object.end)
        self.phase = gvf_feature_line_object.phase # this is always a placeholder'.'

        # VCF DATALINE
        self.chrom = gvf_feature_line_object.seqid
        self.pos = int(gvf_feature_line_object.start)
        self.id = self.vcf_value["ID"]  # attributes: ID
        self.length = self.end - self.pos
        self.qual = gvf_feature_line_object.score # see EVA-3879: this is always '.'
        self.filter = "." # this is always a placeholder '.'; perhaps could add s50.

        # INFO
        #TODO: specific SV info keys populated from gvf_feature_line
        self.key = self.chrom + "_" + str(self.pos)
        self.info = [] # TODO: add info field for self.info
        self.info.append(self.info_string)
        # calculated last
        self.ref = self.get_ref()
        self.alt = self.get_alt(field_lines_dictionary, all_possible_lines_dictionary)

        self.sample_name = self.vcf_value["sample_name"] # this should be each samples names format value # sample names needs to be populated in attributes
        # # higher priority
        self.format = "pending" #TODO: set this in convertgvfattributes
        # # each item in the list exclude_from_info has its own place in the VCF file, so not part of info
        # exclude_from_info = ["ID", # done above
        #                      "Variant_seq", # done above
        #                      "Reference_Seq", # done above
        #                      "Genotype",
        #                      "Phased",
        #                      "sample_name",
        #                      "Zygosity"]
        # for k in self.attributes.keys():
        #     if k not in exclude_from_info:
        #         self.info = self.info + k + "=" + str(self.attributes[k]) + ";"
        #     elif k =="ID":
        #         pass
        #     elif k == "Variant_seq":
        #         #TODO: ADD TO ALT
        #         pass
        #     elif k == "Reference_Seq":
        #         # TODO: ADD TO REF
        #         pass
        #     elif k == "Genotype":
        #         # TODO: ADD TO GT (FORMAT)
        #         pass
        #     elif k == "Phased":
        #         # TODO: ADD TO GT
        #         pass
        #     elif k == "sample_name":
        #         # ADD TO LIST: sample_name
        #         pass
        #     elif k == "Zygosity":
        #         # add to format, gt flag
        #         pass
        #     else:
        #         pass

    def add_padded_base(self, ref, alt, placed_before : bool):
        """ Adds padded base to REF and ALT allele
        :param ref: reference allele
        :param alt: alt allele
        :param placed_before: padded base is placed before ref or alt True or False
        :return: (padded_base, self.pos, self.ref, self.alt)
        """
        if placed_before:
            padded_base_pos = self.pos - 1
            self.pos = padded_base_pos
            padded_base = extract_reference_allele(self.assembly, self.chrom, self.pos, self.end)
            ref = padded_base + ref
            if alt == ".":
                alt = padded_base
            else:
                alt = padded_base + alt
        elif not placed_before:
            padded_base_pos = self.pos + 1
            new_end = self.end + 1
            padded_base = extract_reference_allele(self.assembly, self.chrom, padded_base_pos, new_end)
            ref = ref + padded_base
            if alt == ".":
                alt = padded_base
            else:
                alt = alt + padded_base
        else:
            print("WARNING: Variable placed_before unknown: " + str(placed_before))
            padded_base = None
        return padded_base, self.pos, ref, alt

    def build_iupac_ambiguity_code(self):
        """ Builds dictionary for the iupac ambiguity code
        :return: iupac_ambiguity_dictionary: iupac code as key, list of values as value
        """
        # see PMID: 20202974 (Table 1) for the official list
        iupac_codes = ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]
        R = ["A", "G"]
        Y = ["C", "T"]
        M = ["A", "C"]
        K = ["G", "T"]
        S = ["C", "G"]
        W = ["A", "T"]
        H = ["A", "C", "T"]
        B = ["C", "G", "T"]
        V = ["A", "C", "G"]
        D = ["A", "G", "T"]
        N = ["A", "C", "G", "T"]
        iupac_values = [R, Y, M, K, S, D, W, H, B, V, D, N]
        iupac_ambiguity_dictionary = dict(zip(iupac_codes, iupac_values))
        return iupac_ambiguity_dictionary

    def convert_iupac_ambiguity_code(self, iupac_ambiguity_dictionary, ref_to_convert):
        """ Converts the REF allele if it contains IUPAC ambiguity cod
        :param iupac_ambiguity_dictionary: dictionary of IUPAC ambiguity code and a list of values
        :param ref_to_convert: reference allele to be converted
        :return: self.ref
        """
        converted_ref = ""
        for base in ref_to_convert:
            if base in iupac_ambiguity_dictionary:
                iupac_value = min(iupac_ambiguity_dictionary[base])
                converted_base = iupac_value
            else:
                converted_base = base
            converted_ref = converted_ref + converted_base
        return converted_ref

    def check_ref(self, ref_allele_to_be_checked):
        """ Checks whether a reference allele meets the requirements of the VCF specification
        :param ref_allele_to_be_checked: reference allele to check
        :return: checked_reference_allele: reference allele that meets the requirements of the VCF specification"""
        if isinstance(ref_allele_to_be_checked, str):
            if not all(bases in ref_allele_to_be_checked for bases in ["A", "C", "G", "T", "N"]):
                checked_reference_allele = self.convert_iupac_ambiguity_code(self.iupac_ambiguity_dictionary, ref_allele_to_be_checked)
            else:
                checked_reference_allele = ref_allele_to_be_checked
        else:
            print("WARNING: Ref allele must be a string. Setting to missing value.", ref_allele_to_be_checked)
            checked_reference_allele = "."
        return checked_reference_allele

    def get_ref(self):
        """ Gets the reference allele from attributes column or if not found, returns "."
        :return: reference allele
        """
        if "Reference_seq" in self.vcf_value.keys():
            reference_allele = self.vcf_value["Reference_seq"]
        else:
            if self.assembly:
                reference_allele = extract_reference_allele(self.assembly, self.chrom, self.pos, self.end)
            else:
                print("WARNING: No reference provided. Placeholder inserted for Reference allele.")
                reference_allele = "."
        if reference_allele != ".":
            reference_allele = self.check_ref(reference_allele)
        return reference_allele


    def generate_symbolic_allele(self, field_lines_dictionary, all_possible_lines_dictionary):
        """ Generates the symbolic allele and stores the corresponding metainformation lines. Also determines if variant is precise or imprecise.
        :param field_lines_dictionary: lines for ALT, INFO, etc
        :param all_possible_lines_dictionary: all possible lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        symbolic_allele_id = self.symbolic_allele_dictionary[self.so_type][1]
        symbolic_allele = f'<{symbolic_allele_id}>'

        lines_standard_alt = field_lines_dictionary["ALT"]
        lines_standard_info = field_lines_dictionary["INFO"]
        all_possible_alt_lines = all_possible_lines_dictionary["ALT"]
        all_possible_info_lines = all_possible_lines_dictionary["INFO"]

        if symbolic_allele_id in all_possible_alt_lines:
            lines_standard_alt.append(all_possible_alt_lines[symbolic_allele_id])

        info_svlen = None
        if self.length:
            info_svlen = "SVLEN=" + str(self.length)

        start_range_lower_bound = self.vcf_value["Start_range"][0]
        start_range_upper_bound = self.vcf_value["Start_range"][1]
        end_range_lower_bound = self.vcf_value["End_range"][0]
        end_range_upper_bound = self.vcf_value["End_range"][1]

        # setting up fields to be inserted into INFO
        info_end = None
        info_imprecise = None
        info_cipos = None
        info_ciend = None

        if start_range_lower_bound == "." or start_range_upper_bound == "." or end_range_lower_bound == "." or end_range_upper_bound == ".":
            is_imprecise = False
            info_end = "END=" + str(self.pos + len(self.ref) - 1)
        else:
            is_imprecise = True
            info_imprecise = "IMPRECISE"

            cipos_lower_bound = int(start_range_lower_bound) - self.pos
            cipos_upper_bound = int(start_range_upper_bound) - self.pos
            info_cipos = "CIPOS=" + str(cipos_lower_bound) + "," + str(cipos_upper_bound)

            ciend_lower_bound = int(start_range_lower_bound) - self.pos
            ciend_upper_bound = int(start_range_upper_bound) - self.pos
            info_ciend = "CIEND=" + str(ciend_lower_bound) + "," + str(ciend_upper_bound)

            if symbolic_allele == "<INS>":
                info_end ="END=" + str( self.pos + len(self.ref) - 1 )
            elif symbolic_allele in {"<DEL>", "<DUP>", "<INV>", "<CNV>"}:
                info_end = "END=" + str(self.pos + self.length)
            elif symbolic_allele == "<*>":
                info_end = "END=" + str(self.pos + len(self.ref))
            else:
                "Cannot identify symbolic allele"
        # for all variants (precise and imprecise)
        self.info.append(info_end)
        lines_standard_info.append(all_possible_info_lines["END"])
        self.info.append(info_svlen)
        lines_standard_info.append(all_possible_info_lines["SVLEN"])

        # for imprecise variants only
        if is_imprecise:
            self.info.append(info_imprecise)
            lines_standard_info.append(all_possible_info_lines["IMPRECISE"])
            self.info.append(info_cipos)
            lines_standard_info.append(all_possible_info_lines["CIPOS"])
            self.info.append(info_ciend)
            lines_standard_info.append(all_possible_info_lines["CIEND"])
        return symbolic_allele, self.info, lines_standard_alt, lines_standard_info

    def get_alt(self, field_lines_dictionary, all_possible_lines_dictionary):
        """ Gets the ALT allele for the VCF file
        :param field_lines_dictionary: store INFO,ALT, FILTER, FORMAT lines
        :param all_possible_lines_dictionary: dictionary of all possible ALT, INFO, FORMAT, FILTER lines
        :return: symbolic_allele, self.info, lines_standard_ALT, lines_standard_INFO
        """
        if any(base in self.vcf_value["Variant_seq"] for base in ["A", "C", "G", "T", "N"]):
            alterative_allele = self.vcf_value["Variant_seq"]
        elif self.vcf_value["Variant_seq"] == '.':
            symbolic_allele, self.info, lines_standard_alt, lines_standard_info = self.generate_symbolic_allele(field_lines_dictionary, all_possible_lines_dictionary)
            if symbolic_allele is None:
                alterative_allele = "."
            elif (self.vcf_value["Variant_seq"] == "." or self.vcf_value["Variant_seq"] == "-") and symbolic_allele is not None:
                alterative_allele = symbolic_allele
                # add padded bases
                if self.pos == 1:
                    #print("pos, ref, alt",self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, False)
                    self.ref = self.check_ref(self.ref)
                else:
                    #print("pos, ref, alt", self.pos,self.ref, alterative_allele)
                    padded_base, self.pos, self.ref, self.alt = self.add_padded_base(self.ref, alterative_allele, True)
                    self.ref = self.check_ref(self.ref)
            else:
                alterative_allele = "."
                print("Cannot identify symbolic allele. Variant type is not supported.")
        else:
            alterative_allele = "."
            print("Could not determine the alternative allele.")
        return alterative_allele

    def __str__(self):
        string_to_return = '\t'.join((self.chrom, self.pos, self.key, self.qual, self.filter, self.info, self.source, self.phase, self.end, self.so_type, self.sample_name, self.format))
        return string_to_return
