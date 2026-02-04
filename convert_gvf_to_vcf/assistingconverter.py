"""
This is contains functions to assist the conversion of gvf attributes
"""
import os
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

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
    gvf_attribute_dictionary = get_gvf_attributes(column9_of_gvf)
    vcf_info_values = {} # key is info field value; value is value
    vcf_format_values = {} # key is format field value; value is value
    dropped_gvf_attributes = []
    for attrib_key, attrib_value in gvf_attribute_dictionary.items():
        if attrib_key in mapping_attribute_dict:
            field_values = mapping_attribute_dict[attrib_key]
            for field in field_values:
                if field == "INFO":
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
                elif field == "FORMAT":
                    field_lines_dictionary[field].append(all_possible_lines_dictionary[field][field_values[field]["FieldKey"]])
                    sample_name = gvf_attribute_dictionary.get("sample_name")
                    if sample_name in vcf_format_values:
                        vcf_format_values[sample_name].update({field_values[field]["FieldKey"]: gvf_attribute_dictionary[attrib_key]})
                    else:
                        vcf_format_values[sample_name] = {field_values[field]["FieldKey"]: gvf_attribute_dictionary[attrib_key]}
                else:
                    logger.warning(f"Unsupported Field: {field}")
        else:
            logger.info(f"catching attribute keys for review at a later date {attrib_key} {attrib_value}")
            dropped_gvf_attributes.append(attrib_key)
    return gvf_attribute_dictionary, vcf_info_values, vcf_format_values
