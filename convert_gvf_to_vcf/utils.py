"""This contains readers and utilities"""
import importlib
import json
import sys
from importlib.metadata import distributions

import yaml
import logging

from convert_gvf_to_vcf.gvf_feature import GvfFeatureline

logger = logging.getLogger(__name__)


def read_in_json_schema(json_schema_file):
    """Reads in JSON file and returns JSON schema as a dictionary
    :param json_schema_file: path to schema
    :return: schema: dictionary of JSON schema
    """
    with open(json_schema_file) as json_file:
        schema = json.load(json_file)
    return schema

def read_yaml(yaml_file):
    """Reads a yaml file (of attributes) and returns dictionary
    :param yaml_file: file of attributes
    :return: : A dictionary of attributes.
    """
    with open(yaml_file) as y_file:
        attributes = yaml.safe_load(y_file)
    return attributes

def read_pragma_mapper(pragma_mapper_file):
    """ Reads in the pragma mapper file and stores as dictionary.
    :param pragma_mapper_file: A file containing the pragma and their VCF equivalent
    :return: pragma_to_vcf_header: dictionary of pragma => vcf header
    """
    pragma_to_vcf_header = {}
    with open(pragma_mapper_file) as pragma_file:
        next(pragma_file)
        for line in pragma_file:
            pragma_tokens = line.rstrip().split("\t")
            pragma = pragma_tokens[0]
            vcf_header = pragma_tokens[1]
            pragma_to_vcf_header[pragma] = vcf_header
    return pragma_to_vcf_header

def read_in_gvf_header(gvf_input):
    """ Reads in the user provided GVF file.
    :param gvf_input: The input GVF file
    :return:
        - gvf_pragmas: list of pragma lines (start with ## at the top of GVF file)
        - gvf_non_essential: list of non essential pragma (start with # near the top of GVF file)
    """
    gvf_pragmas = []  # list of pragma lines starting with: ##
    gvf_non_essential = []  # list of non-essential lines starting with: #
    with open(gvf_input) as gvf_file:
        for line in gvf_file:
            if line.startswith("##"):
                gvf_pragmas.append(line.rstrip())
            elif line.startswith("#"):
                gvf_non_essential.append(line.rstrip())
            else:
                break
    return gvf_pragmas, gvf_non_essential

def read_in_gvf_data(gvf_input):
    """ Reads in the user provided GVF file.
    :param gvf_input: arguments.gvf_input : The input GVF file
    :return iterator of GvfFeatureline objects
    """

    with open(gvf_input) as gvf_file:
        for line in gvf_file:
            if not line.startswith("#"):
                f_list = line.rstrip().split("\t")
                yield GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])

def generate_symbolic_allele_dict(mapping_dictionary):
    """Reads in mapping dictionary and returns a symbolic allele dictionary.
    :param mapping_dictionary: mapping dictionary
    :return symbolic_allele_dict: stores information for a symbolic allele
    """
    symbolic_allele_dict = {}
    for attribute in mapping_dictionary:
        # Symbolic alleles refer only to the header type "ALT".
        header_type = "ALT"
        if mapping_dictionary[attribute].get(header_type) is not None:
            if mapping_dictionary[attribute].get(header_type).get("FieldKey") is not None:
                name = attribute
                sequence_ontology_id = mapping_dictionary[attribute].get(header_type).get("SOID")
                symb_allele = mapping_dictionary[attribute].get(header_type).get("FieldKey")
                description = mapping_dictionary[attribute].get(header_type).get("Description")

                symbolic_allele_dict.setdefault(name, []).append(sequence_ontology_id)
                symbolic_allele_dict.setdefault(name, []).append(symb_allele)
                symbolic_allele_dict.setdefault(name, []).append(description)
    return symbolic_allele_dict

def build_iupac_ambiguity_code():
    """ Builds dictionary for the iupac ambiguity code.
    :return: iupac_ambiguity_dictionary: iupac code as key, list of values as value
    """
    # see PMID: 20202974 (Table 1) for the official list
    return {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "M": ["A", "C"],
        "K": ["G", "T"],
        "S": ["C", "G"],
        "W": ["A", "T"],
        "H": ["A", "C", "T"],
        "B": ["C", "G", "T"],
        "V": ["A", "C", "G"],
        "D": ["A", "G", "T"],
        "N": ["A", "C", "G", "T"]
    }

def get_validated_value(config, key_parts_to_get, expected_type, default_value=None):
    """Gets information from a config file (ebi_eva_common_pyutils.config > cfg.load_config_file)
    :param config: config
    :param key_parts_to_get: key to retrieve from config file (hierarchical)
    :param expected_type: expected_type of key to retrieve from config file
    :param default_value: default_value of key to retrieve from config file
    :return: value_in_config
    """
    # key_parts_to_get is a tuple of multiple values so unpacking
    value_in_config = config.query(*key_parts_to_get, ret_default=default_value)
    if value_in_config is None:
        raise ValueError(f"Missing key required: {key_parts_to_get}")
    if not isinstance(value_in_config, expected_type):
        try:
            # cast the value to its expected type
            value_in_config = expected_type(value_in_config)
        except (ValueError, TypeError):
            raise TypeError(f"Key '{key_parts_to_get}' must be {expected_type.__name__}")
    return value_in_config

def get_version():
    """Gets the dynamic version of the name in pyproject.toml file"""
    try:
        return importlib.metadata.version("convertGVFtoVCF")
    except Exception:
        return "Not found"

def log_environment_packages(logger, file_name, VERSION):
    """Logs the versions of packages used."""
    logger.info(f"===={file_name} v{VERSION}====")
    logger.info(f"Running on Python version: {sys.version}")
    logger.info("==== FULL INSTALLED PACKAGES LIST ====")

    # Extracts every single library installed in this environment
    all_packages = sorted({f"{d.metadata['Name']}=={d.version}" for d in distributions()})

    for package in all_packages:
        logger.info(f"[PACKAGE] {package}")
    logger.info("=======================================")