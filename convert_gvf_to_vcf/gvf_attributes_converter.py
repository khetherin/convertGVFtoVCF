"""
This contains functions to assist the conversion of gvf attributes
"""
from ebi_eva_common_pyutils.logger import logging_config as log_cfg


logger = log_cfg.get_logger(__name__)

def generate_custom_structured_meta_line(field, idkey, number,
                                         data_type, description,
                                         optional_data=None):
    """ Generates a custom structured meta-information line for INFO/FILTER/FORMAT/ALT
    :param field: required field INFO, FILTER, FORMAT, ALT
    :param idkey: required field for structured lines ID
    :param number: Number of values included or special character: A or R or G or .
    :param data_type: Values are Integer, Float, Character, String
    :param description: Description
    :param optional_data: an optional field, dictionary of custom fields and their values
    :return: custom_structured_string
    """
    extra_keys_kv_lines = []
    if optional_data:
        for extra_field in optional_data:
            kv_line = "," + extra_field + "=" + '"' + optional_data[extra_field] + '"'
            extra_keys_kv_lines.append(kv_line)
    vcf_extra_keys = ''.join(extra_keys_kv_lines)
    custom_structured_string = (f'##{field}=<'
                                f'ID={idkey},'
                                f'Number={number},'
                                f'Type={data_type},'
                                f'Description="{description}"'
                                f'{vcf_extra_keys}>')
    return custom_structured_string

class GvfAttributeParser:
    """The responsibility of this class is to parse the attributes in column 9 of the GVF file"""
    def __init__(self, column9_of_gvf):
        """Initialise the parser and extract gvf attributes
        :param column9_of_gvf: column of the GVF file
        """
        self._raw_column9 = column9_of_gvf
        self._attributes = self.get_gvf_attributes(column9_of_gvf)

    def get_gvf_attributes(self, column9_of_gvf):
        """Get a dictionary of GVF attributes
        :param column9_of_gvf:  column - the final column of the GVF file
        :return: gvf_attribute_dictionary: a dictionary of attribute keys and their values
        """
        gvf_attribute_dictionary = {}  # attribute key => value
        # parse by semicolon this creates attribute
        # parse by equals sign this creates tag-values, if the value is a comma, create a list
        attributes_in_gvf_line = column9_of_gvf.split(";")
        for attribute in attributes_in_gvf_line:
            attribute_key, attribute_value = attribute.split("=")
            if "," in attribute_value:
                attribute_value_list = attribute_value.split(",")
                gvf_attribute_dictionary[attribute_key] = attribute_value_list
            else:
                gvf_attribute_dictionary[attribute_key] = attribute_value
        return gvf_attribute_dictionary

class GvfAttributeTransformer:
    """The responsibility of this class is to convert the GVF attributes to a VCF-compliant fields."""
    def __init__(self, mapping_attribute_dict, field_lines_dictionary, all_possible_lines_dictionary):
        """Initialises the transformer with the dictionaries it  needs to convert GVF attributes
        :param mapping_attribute_dict: Dictionary mapping GVF attributes to VCF fields.
        :param field_lines_dictionary: Dictionaries for tracking header lines (ALT, INFO, FILTER, and FORMAT).
        :param all_possible_lines_dictionary: All possible VCF header lines.
        """
        self._mapping_attribute_dict = mapping_attribute_dict
        self._field_lines_dictionary = field_lines_dictionary
        self._all_possible_lines_dictionary = all_possible_lines_dictionary

    def convert_gvf_attributes_to_vcf_values(self,
                                             gvf_parser):
        """Converts GVF attributes to a dictionary that will store VCF values.
        Populates ALT INFO FILTER FORMAT with the correct VCF values.
        :param gvf_parser: GvfAttributeParser object to obtain attributes column 9 of gvf file
        :return gvf_attribute_dictionary, vcf_info_values, vcf_format_values: dict of GVF attributes and VCF values.
        """
        # parse GVF attributes
        gvf_attribute_dictionary = gvf_parser._attributes
        # create dictionaries to store INFO and FORMAT values
        vcf_info_values = {} # key is info field value; value is value
        vcf_format_values = {} # key is format field value; value is value
        dropped_gvf_attributes = [] # create list to track unmapped attributes that do not belong in mapping_attribute_dict
        for attrib_key, attrib_value in gvf_attribute_dictionary.items():
            # determine if this found in the mapping rules
            if attrib_key in self._mapping_attribute_dict:
                field_values = self._mapping_attribute_dict[attrib_key]
                for field in field_values:
                    self.process_vcf_fields(attrib_key, field, field_values, gvf_parser,
                                            vcf_format_values, vcf_info_values)
            else:
                logger.info(f"catching attribute keys for review at a later date {attrib_key} {attrib_value}")
                dropped_gvf_attributes.append(attrib_key)
        self.apply_genotype_inference(gvf_parser, vcf_format_values, vcf_info_values)
        return gvf_attribute_dictionary, vcf_info_values, vcf_format_values

    def apply_genotype_inference(self, gvf_parser, vcf_format_values, vcf_info_values):
        """Determines and infers genotype values if the VCF record is SITES-ONLY.
        :param gvf_parser: GvfAttributeParser object
        :param vcf_format_values: dictionary of sample names to FORMAT field key-value pairs
        :param vcf_info_values: dictionary of sample names to INFO field key-value pairs
        """
        vcf_type = self.determine_vcf_type(vcf_info_values, vcf_format_values)
        if vcf_type == "SITES-ONLY":
            self.infer_genotype(gvf_parser._attributes, vcf_format_values)

    def process_vcf_fields(self, attrib_key, field, field_values, gvf_parser, vcf_format_values,
                           vcf_info_values):
        """Process depending on VCF field.
        :param attrib_key: gvf attribute to be processed: key to extract value from gvf_attribute_dictionary
        :param field: current field to be processed
        :param field_values: configuration dictionary (FieldKey, Number, Type, Description) by field
        :param gvf_parser: GvfAttributeParser object
        :param vcf_format_values: dictionary of sample names to FORMAT field key-value pairs
        :param vcf_info_values: dictionary of sample names to INFO field key-value pairs
        """
        if field == "INFO":
            self.process_vcf_info_field(attrib_key, field, field_values, gvf_parser._attributes, vcf_info_values)
        elif field == "FORMAT":
            self.process_vcf_format_field(attrib_key, field, field_values, gvf_parser._attributes, vcf_format_values)
        else:
            logger.warning(f"Unsupported Field: {field}")


    def infer_genotype(self, gvf_attribute_dictionary, vcf_format_values):
        sample_name = gvf_attribute_dictionary.get("sample_name") or "UNKNOWN_SAMPLE"
        if sample_name not in vcf_format_values:
            vcf_format_values[sample_name] = {}
        if "GT" not in vcf_format_values[sample_name]:
            vcf_format_values[sample_name]["GT"] = "1/."


    def determine_vcf_type(self, vcf_info_values, vcf_format_values):
        """Determines the type of VCF for the VCF line.
        Uses INFO and FORMAT values to determine if genotyped, populated or site-only.
        :param vcf_info_values: dictionary of INFO values
        :param vcf_format_values: dictionary of FORMAT values
        """
        has_genotypes = any(isinstance(sample_data, dict) and "GT" in sample_data
                            for sample_data in vcf_format_values.values())
        has_allele_frequencies = "AF" in vcf_info_values
        has_allele_counts = ("AC" in vcf_info_values) and ("AN" in vcf_info_values)
        if has_genotypes:
            return "GENOTYPED"
        elif has_allele_frequencies and has_allele_counts:
            return "POPULATED"
        # if not an EVA-accepted VCF file (see https://www.ebi.ac.uk/eva/?Submit-Data)
        elif not ( has_genotypes or has_allele_frequencies or has_allele_counts):
            logger.info("The VCF file has been determined as site-only. This will be converted to genotypes. ")
            return "SITES-ONLY"
        else:
            logger.info("unable to infer the vcf type")
            return "UNKNOWN"

    def process_vcf_format_field(self, attrib_key, field, field_values,
                                 gvf_attribute_dictionary, vcf_format_values):
        """Maps GVF attributes to VCF FORMAT fields and stores them in a dictionary
        :param attrib_key: key to extract value from gvf_attribute_dictionary
        :param field: current field to be processed
        :param field_values: configuration dictionary (FieldKey, Number, Type, Description) by field
        :param gvf_attribute_dictionary: parsed attributes from the GVF line
        :param vcf_format_values: dictionary of sample names to FORMAT field key-value pairs
        """
        target_key = field_values[field]["FieldKey"]
        self._field_lines_dictionary[field].append(self._all_possible_lines_dictionary[field][target_key])
        sample_name = gvf_attribute_dictionary.get("sample_name")
        if sample_name in vcf_format_values:
            vcf_format_values[sample_name].update({target_key: gvf_attribute_dictionary[attrib_key]})
        else:
            vcf_format_values[sample_name] = {target_key: gvf_attribute_dictionary[attrib_key]}


    def process_vcf_info_field(self, attrib_key, field, field_values, gvf_attribute_dictionary,
                               vcf_info_values):
        """Maps GVF attributes to VCF INFO field and stores them in a dictionary
        :param attrib_key: key to extract value from gvf_attribute_dictionary
        :param field: current field to be processed
        :param field_values: configuration dictionary (FieldKey, Number, Type, Description) by field
        :param gvf_attribute_dictionary: parsed attributes from the GVF line
        :param vcf_info_values: dictionary of sample names to INFO field key-value pairs
        """
        header = generate_custom_structured_meta_line(
            field=field,
            idkey=field_values[field]["FieldKey"],
            number=field_values[field]["Number"],
            data_type=field_values[field]["Type"],
            description=field_values[field]["Description"],
            optional_data=None
        )
        self._field_lines_dictionary[field].append(header)
        vcf_info_values[field_values[field]["FieldKey"]] = gvf_attribute_dictionary[attrib_key]


