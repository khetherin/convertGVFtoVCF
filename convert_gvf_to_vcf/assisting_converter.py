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

def get_gvf_attributes(column9_of_gvf):
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


def convert_gvf_attributes_to_vcf_values(column9_of_gvf,
                                         mapping_attribute_dict,
                                         field_lines_dictionary,
                                         all_possible_lines_dictionary):
    """Converts GVF attributes to a dictionary that will store VCF values.
    Populates ALT INFO FILTER FORMAT with the correct VCF values.
    :param column9_of_gvf: attributes column of gvf file
    :param mapping_attribute_dict: attributes
    :param field_lines_dictionary: dictionaries for ALT INFO FILTER and FORMAT
    :param all_possible_lines_dictionary: all possible VCF header lines
    :return gvf_attribute_dictionary, vcf_info_values, vcf_format_values: dict of GVF attributes and VCF values.
    """
    # parse GVF attributes
    gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
    # create dictionaries to store INFO and FORMAT values
    vcf_info_values = {} # key is info field value; value is value
    vcf_format_values = {} # key is format field value; value is value
    # create list to track unmapped attributes that do not belong in mapping_attribute_dict
    dropped_gvf_attributes = []
    for attrib_key, attrib_value in gvf_attribute_dictionary.items():
        # determine if this found in the mapping rules
        if attrib_key in mapping_attribute_dict:
            field_values = mapping_attribute_dict[attrib_key]
            # determine if this is found in the INFO or FORMAT COLUMN
            for field in field_values:
                if field == "INFO":
                    process_vcf_info_field(attrib_key, field, field_lines_dictionary, field_values,
                                           gvf_attribute_dictionary, vcf_info_values)
                elif field == "FORMAT":
                    process_vcf_format_field(all_possible_lines_dictionary, attrib_key, field, field_lines_dictionary,
                                             field_values, gvf_attribute_dictionary, vcf_format_values)
                else:
                    logger.warning(f"Unsupported Field: {field}")
        else:
            logger.info(f"catching attribute keys for review at a later date {attrib_key} {attrib_value}")
            dropped_gvf_attributes.append(attrib_key)
    vcf_type = determine_vcf_type(vcf_info_values,vcf_format_values)
    if vcf_type == "SITES-ONLY":
        infer_genotype(gvf_attribute_dictionary, vcf_format_values)

    return gvf_attribute_dictionary, vcf_info_values, vcf_format_values

def infer_genotype(gvf_attribute_dictionary, vcf_format_values):
    sample_name = gvf_attribute_dictionary.get("sample_name") or "UNKNOWN_SAMPLE"
    if sample_name not in vcf_format_values:
        vcf_format_values[sample_name] = {}
    if "GT" not in vcf_format_values[sample_name]:
        vcf_format_values[sample_name]["GT"] = "1/."


def determine_vcf_type(vcf_info_values, vcf_format_values):
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

def process_vcf_format_field(all_possible_lines_dictionary, attrib_key, field, field_lines_dictionary, field_values,
                             gvf_attribute_dictionary, vcf_format_values):
    """Maps GVF attributes to VCF FORMAT fields and stores them in a dictionary
    :param all_possible_lines_dictionary: a master dictionary of all possible lines
    :param attrib_key: key to extract value from gvf_attribute_dictionary
    :param field: current field to be processed
    :param field_lines_dictionary: dictionary of compiled VCF lines per field (value is list)
    :param field_values: configuration dictionary (FieldKey, Number, Type, Description) by field
    :param gvf_attribute_dictionary: parsed attributes from the GVF line
    :param vcf_format_values: dictionary of sample names to FORMAT field key-value pairs
    """
    field_lines_dictionary[field].append(all_possible_lines_dictionary[field][field_values[field]["FieldKey"]])
    sample_name = gvf_attribute_dictionary.get("sample_name")
    if sample_name in vcf_format_values:
        vcf_format_values[sample_name].update({field_values[field]["FieldKey"]: gvf_attribute_dictionary[attrib_key]})
    else:
        vcf_format_values[sample_name] = {field_values[field]["FieldKey"]: gvf_attribute_dictionary[attrib_key]}


def process_vcf_info_field(attrib_key, field, field_lines_dictionary, field_values, gvf_attribute_dictionary,
                           vcf_info_values):
    """Maps GVF attributes to VCF INFO field and stores them in a dictionary
    :param attrib_key: key to extract value from gvf_attribute_dictionary
    :param field: current field to be processed
    :param field_lines_dictionary: dictionary of compiled VCF lines per field (value is list)
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
    field_lines_dictionary[field].append(header)
    vcf_info_values[field_values[field]["FieldKey"]] = gvf_attribute_dictionary[attrib_key]
