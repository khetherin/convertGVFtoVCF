import argparse
import os
from collections import Counter
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from convert_gvf_to_vcf.conversionstatistics import StatisticsPayload
from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.utils import read_in_gvf_header, read_in_gvf_data
from convert_gvf_to_vcf.vcfline import VcfLineBuilder


logger = log_cfg.get_logger(__name__)

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

# the functions below relate to the VCF header (Part 1)
def generate_vcf_header_structured_lines(header_type, mapping_attribute_dict):
    """ Generates a dictionary of all possible standard structured lines for INFO/FILTER/FORMAT/ALT.
    :param header_type: type of header file to read i.e. ALT, FILTER, INFO or FORMAT
    :param mapping_attribute_dict: dictionary of all attributes
    :return: dictionary of all possible standard structured lines keys for the header type
    """
    all_possible_lines = {}

    for attribute in mapping_attribute_dict:
        # Formatting the header string for FILTER, INFO or FORMAT and storing in a dictionary
        if mapping_attribute_dict[attribute].get(header_type) is not None and header_type != "ALT":
            header_string = (f'##{header_type}='
                             f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                             f'Number={mapping_attribute_dict[attribute][header_type]["Number"]},'
                             f'Type={mapping_attribute_dict[attribute][header_type]["Type"]},'
                             f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
            all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        # Formatting the header string for ALT and storing in a dictionary
        elif mapping_attribute_dict[attribute].get(header_type) is not None and header_type == "ALT":
            if mapping_attribute_dict[attribute][header_type]["FieldKey"] is not None:
                header_string = (f'##{header_type}='
                                 f'<ID={mapping_attribute_dict[attribute][header_type]["FieldKey"]},'
                                 f'Description="{mapping_attribute_dict[attribute][header_type]["Description"]}">')
                all_possible_lines[mapping_attribute_dict[attribute][header_type]["FieldKey"]] = header_string
        else:
            pass
    return all_possible_lines

def generate_vcf_header_unstructured_line(vcf_unstructured_key,
                                          vcf_unstructured_value):
    """ Generates a formatted unstructured metainformation line using a custom key value pair e.g. "##key=value"
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    return custom_unstructured_string

def get_sample_name_from_pragma(pragma_name, pragma_value):
    """Returns the sample name in a list from a pragma name and value
    :param pragma_name: First part of the pragma e.g.: #sample
    :param pragma_value: Second part of the pragma e.g.: sample_name=Wilds2-3;subject_name=Wilds2-3
    :return: sample_names: a list (expect to contain one element - to be extended later on )
    """
    sample_name = None
    if pragma_name.startswith("#sample"):
        list_of_sample_information = pragma_value.split(";")
        for sample_info in list_of_sample_information:
            if sample_info.startswith("sample_name"):
                sample_name = sample_info.split("=")[1]
    return sample_name

def get_unique_sample_names(sample_names):
    """ Takes a list of sample names and returns a unique list
    :param: sample_names: list of sample names (can contain duplicates)
    :return: uniq_sample_name: list of unique sample names
    """
    seen_sample_names = set()
    uniq_sample_name = []
    for sample in sample_names:
        if sample not in seen_sample_names:
            seen_sample_names.add(sample)
            uniq_sample_name.append(sample)
    return uniq_sample_name


def convert_gvf_pragmas_to_vcf_header(list_of_gvf_pragmas_to_convert,
                                      list_of_gvf_pragmas,
                                      pragma_to_vcf_map):
    """ This function converts pragmas using the delimiter for pragmas (delim=" ").
    Format of pragma = Pragmas are found in the header and must start with "##".
    Function of pragmas = Pragmas are essential for compliance with the specification and to be considered a valid GVF file. Important for interoperability.
    :param: list_of_gvf_pragmas_to_convert: pragmas to be converted
    :param: list_of_converted_pragmas: list of existing converted pragmas: will be appending results to this list
    :param: list_of_essential_pragma: a reference list of essential pragmas
    :param: pragma_to_vcf_map: a mapping dict of GVF pragmas and their VCF counterpart
    : return: list_of_converted_pragmas: updated version
    """
    list_of_converted_pragmas = []
    for pragma in list_of_gvf_pragmas_to_convert:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma, " ", list_of_gvf_pragmas, pragma_to_vcf_map)
        list_of_converted_pragmas.append(generate_vcf_header_unstructured_line(vcf_header_key, pragma_value))
    return list_of_converted_pragmas

def convert_gvf_pragma_comment_to_vcf_header(gvf_pragma_comments_to_convert,
                                             list_of_gvf_pragma_comments,
                                             pragma_to_vcf_map):
    """ This converts pragma comments (nonessential) including obtaining the sample names from the gvf pragma
    Format of pragma comment = These tend to start with '#'. These comments are relevant to DGVa delimiter is ": "
    Function of pragma comment = These are non-essential and tend to be ignored by GVF processors. They can contain useful additional info.
    :param: gvf_pragma_comments_to_convert: pragma comments to be converted
    :param: list_of_converted_pragma_comments: will append results to this list
    :param: list_of_gvf_pragma_comments : reference list
    :param: pragma_to_vcf_map: a mapping dict of GVF pragmas and their VCF counterpart
    :param: sample_names_from_pragma_comments: will append results to this list
    return: list_of_converted_pragma_comments, sample_names_from_pragma_comments
    """
    list_of_converted_pragma_comments = []
    sample_names_from_pragma_comments = []
    for gvf_pragma_comment in gvf_pragma_comments_to_convert:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(gvf_pragma_comment, ": ", list_of_gvf_pragma_comments, pragma_to_vcf_map)
        if pragma_name.startswith("#Publication"):
            if ";" in pragma_value:
                publication_tokens = get_pragma_tokens(pragma_value, ";", "=")
                for pub_token in publication_tokens:
                    list_of_converted_pragma_comments.append(generate_vcf_header_unstructured_line(pub_token[0], pub_token[1]))
            else:
                list_of_converted_pragma_comments.append(generate_vcf_header_unstructured_line(pragma_name.lstrip("#"), pragma_value))
        elif pragma_name == "#Study":
            study_tokens = get_pragma_tokens(pragma_value, ";", "=")
            for s_token in study_tokens:
                list_of_converted_pragma_comments.append(generate_vcf_header_unstructured_line(s_token[0], s_token[1]))
        else:
            if vcf_header_key is not None:
                list_of_converted_pragma_comments.append(generate_vcf_header_unstructured_line(vcf_header_key, pragma_value))
        # populating sample headers
        sample_name = get_sample_name_from_pragma(pragma_name, pragma_value)
        if sample_name is not None:
            sample_names_from_pragma_comments.append(sample_name)
    return list_of_converted_pragma_comments, sample_names_from_pragma_comments


def convert_gvf_pragmas_for_vcf_header(gvf_pragmas,
                                       gvf_pragma_comments,
                                       reference_lookup):
    """ Generates a list of VCF-compatible pragmas and samples.
    :param gvf_pragmas: list of gvf pragmas to convert
    :param gvf_pragma_comments: list of gvf pragma comments to convert
    :return: unique_converted_pragmas, sample_names: a list of pragmas (removed duplicates), list of sample names
    """
    converted_pragmas = []
    unique_converted_pragmas = []
    # sample_names = []
    ####
    # MANDATORY header for VCF: file format for VCF
    mandatory_first_line = generate_vcf_header_unstructured_line("fileformat", "VCFv4.4")
    #Go through pragmas
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    list_of_gvf_pragmas = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
    list_of_converted_pragmas = convert_gvf_pragmas_to_vcf_header(gvf_pragmas, list_of_gvf_pragmas, reference_lookup.pragma_to_vcf_map)
    converted_pragmas.extend(list_of_converted_pragmas)
    # Go through gvf pragma comments
    list_of_gvf_pragma_comments = ["#sample", "#Study_accession", "#Study_type", "#Display_name", "#Publication",
                                    "#Study", "#Assembly_name", "#subject"]
    list_of_converted_pragma_comments, sample_names = convert_gvf_pragma_comment_to_vcf_header(gvf_pragma_comments, list_of_gvf_pragma_comments, reference_lookup.pragma_to_vcf_map)
    converted_pragmas.extend(list_of_converted_pragma_comments)

    # ensure unique sample names and preserve order
    unique_sample_name = get_unique_sample_names(sample_names)
    ###
    # unique_converted_pragmas = list(dict.fromkeys(pragma for pragma in converted_pragmas if pragma not in unique_converted_pragmas))
    unique_converted_pragmas = list(set(converted_pragmas))
    unique_converted_pragmas = [mandatory_first_line] + unique_converted_pragmas
    return unique_converted_pragmas, unique_sample_name

# the function below relates to the VCF headerline (Part 2)
def generate_vcf_header_line(is_missing_format, samples):
    """ Generates the VCF header line using the nine mandatory headers and the sample names.
    :param is_missing_format: boolean (False = Format + sample names, True = mandatory fields only)
    :param samples: list of samples, these will appear in the header line
    :return: vcf_header: a string
    """
    if is_missing_format:
        vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        vcf_header = '\t'.join(vcf_header_fields)
    else:

        vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        for sample in samples:
            vcf_header_fields.append(sample)
        vcf_header = '\t'.join(vcf_header_fields)
    return vcf_header

# the functions below relate to the GVF header
def parse_pragma(pragma_to_parse, delimiter):
    """ Parses pragma and returns name and value of the pragma.
    :param pragma_to_parse: pragma
    :param delimiter: to split by
    :return: pragma_name, pragma_value: key and value of pragma
    """
    try:
        pragma_tokens = pragma_to_parse.split(delimiter)
        pragma_name = pragma_tokens[0]
        if len(pragma_tokens) >= 2:
            pragma_value = ''.join(map(str, pragma_tokens[1:]))
        # elif len(pragma_tokens) == 1:
            # pragma_value = ''.join(map(str, pragma_tokens[0]))
        else:
            pragma_value = None
            logger.warning(f"WARNING: no value for the following pragma {pragma_to_parse}")
        return pragma_name, pragma_value
    except AttributeError as e:
        logger.error(f"Skipping this, can't be parsed {pragma_to_parse}: {e}")
        raise AttributeError(f"Cannot parse {pragma_to_parse}")

def get_pragma_name_and_value(pragma_to_parse, delimiter, pragma_list, pragma_name_to_vcf_dict):
    """Get pragma name and value and its corresponding VCF header key.
    :param pragma_to_parse: pragma that will be parsed
    :param delimiter: the separator
    :param pragma_list: list of pragmas to search through
    :param pragma_name_to_vcf_dict: dictionary pragma name and its vcf entry
    :return vcf_header_key, pragma_name, pragma_value
    """
    pragma_name, pragma_value = parse_pragma(pragma_to_parse, delimiter)
    if pragma_name in pragma_list:
        vcf_header_key = pragma_name_to_vcf_dict.get(pragma_name)
    else:
        vcf_header_key = None
    return vcf_header_key, pragma_name, pragma_value

def get_pragma_tokens(pragma_value, first_delimiter, second_delimiter):
    """Get pragma tokens for nested pragmas
    :param pragma_value: value to parse
    :param first_delimiter: first separator
    :param second_delimiter: second separtor
    :return pragma_tokens
    """
    initial_list = pragma_value.split(first_delimiter)
    pragma_tokens = []
    for element in initial_list:
        pragma_tokens.append(element.split(second_delimiter))
    return pragma_tokens


def write_header(vcf_output, pragmas_for_vcf, header_lines_per_type, header_fields, samples):
    logger.info(f"Total number of samples in this VCF: {len(samples)}")
    vcf_header_file = vcf_output + '_header'
    with open(vcf_header_file, "w") as vcf_header_output:

        # Part 1 of VCF file: Write the VCF header. This will include preserved data from the GVF file.
        for pragma in pragmas_for_vcf:
            vcf_header_output.write(f"{pragma}\n")
        for header_type in header_lines_per_type.values():
            for header_line in header_type:
                vcf_header_output.write(f"{header_line}\n")

        # Part 2 of VCF file: Write the VCF header line.
        # Write the header.
        vcf_header_output.write(f"{header_fields}\n")
    return vcf_header_file

def convert(gvf_input, vcf_output, assembly):
    # Log the inputs and outputs.
    logger.info("Running the GVF to VCF converter")
    logger.info(f"The provided input file is: {gvf_input}")
    logger.info(f"The provided output file is: {vcf_output}")
    if assembly:
        logger.info(f"The provided assembly file is: {assembly}")
    assembly_file = os.path.abspath(assembly)
    assert os.path.isfile(assembly_file), "Assembly file does not exist"

    # Read input file and separate out its components
    logger.info(f"Reading in the following GVF header from {gvf_input}")
    gvf_pragmas, gvf_pragma_comments = read_in_gvf_header(gvf_input)
    # Creating lookup object to store important dictionaries and log what has been stored.
    reference_lookup = Lookup(assembly_file)
    logger.info("Creating the reference lookup object.")
    logger.info("Storing the attributes file: attribute_mapper.yaml")
    logger.info("Storing the symbolic allele dictionary.")
    logger.info(f"Storing the assembly file: {assembly_file}")
    logger.info("Storing the IUPAC ambiguity dictionary.")

    # Preparation work:
    # Store the VCF metainformation and ensure preservation of important GVF data.
    # This information will be useful when creating the VCF header.
    # TODO: refactor function generate_vcf_metainfo
    (
        pragmas_for_vcf,
        samples
    ) = convert_gvf_pragmas_for_vcf_header(gvf_pragmas, gvf_pragma_comments, reference_lookup)

    # TODO: place the all_header_lines_per_type_dict into the reference_lookup.

    # Create data structure to store all possible outcomes for header lines (for fields ALT, INFO, FILTER, FORMAT)
    all_header_lines_per_type_dict = {
        htype: generate_vcf_header_structured_lines(htype, reference_lookup.mapping_attribute_dict) for htype in
        ["ALT", "INFO", "FILTER", "FORMAT"]
    }
    vcf_builder = VcfLineBuilder(all_header_lines_per_type_dict, reference_lookup, samples)
    is_missing_format_value = True
    # Convert each feature line in the GVF file to a VCF object (stores all the data for a line in the VCF file).
    # NOTE: Main Logic lives here.
    gvf_chromosome_counter = Counter()
    vcf_data_file = vcf_output + '_data_lines'
    with open(vcf_data_file, "w") as open_data_lines:
        logger.info("Generating the VCF datalines")
        previous_vcf_line = None
        gvf_feature_line_count = 0
        for gvf_lines_obj in read_in_gvf_data(gvf_input):
            gvf_feature_line_count += 1
            gvf_chromosome_counter.update([gvf_lines_obj.seqid])
            current_vcf_line = vcf_builder.build_vcf_line(gvf_lines_obj)
            # is_missing_format_value will only be true if all the format field are missing.
            is_missing_format_value = is_missing_format_value and current_vcf_line.format_keys == ['.']
            # Each GVF feature has been converted to a VCF object so begin comparing and merging the VCF objects.
            if previous_vcf_line:
                if current_vcf_line == previous_vcf_line:
                    current_vcf_line.merge(previous_vcf_line, list_of_sample_names=samples)
                else:
                    open_data_lines.write(str(previous_vcf_line) + "\n")
            previous_vcf_line = current_vcf_line
        if previous_vcf_line:
            open_data_lines.write(str(previous_vcf_line) + "\n")
        else:
            logger.warning("No feature lines were found for this GVF file.")

    header_lines_per_type = vcf_builder.build_vcf_header()
    header_fields = generate_vcf_header_line(is_missing_format=is_missing_format_value, samples=samples)
    # vcf_header_file = write_header(vcf_output, pragmas_for_vcf, header_lines_per_type,
                                   # is_missing_format_value, samples)
    vcf_header_file = write_header(vcf_output, pragmas_for_vcf, header_lines_per_type,
                                   header_fields, samples)

    logger.info(f"Combining the header and data lines to the following VCF output: {vcf_output}")
    vcf_chromosome_counter = Counter()
    vcf_data_line_count = 0
    with open(vcf_output, "w") as vcf_output:
        with open(vcf_header_file, "r") as vcf_header_fh:
            for line in vcf_header_fh:
                vcf_output.write(line)
        with open(vcf_data_file, "r") as vcf_data_fh:
            for line in vcf_data_fh:
                vcf_output.write(line)
                vcf_data_line_count += 1
                vcf_chromosome_counter.update([line.split("\t")[0]])
    vcf_output.close()
    logger.info("Remove the temporary files")
    if os.path.exists(vcf_header_file):
        os.remove(vcf_header_file)
    if os.path.exists(vcf_data_file):
        os.remove(vcf_data_file)
    logger.info("GVF to VCF conversion complete")
    conversion_payload = StatisticsPayload(
        samples=samples,
        vcf_header_fields=header_fields,
        gvf_line_count= gvf_feature_line_count,
        vcf_line_count= vcf_data_line_count,
        gvf_chromosome_count=dict(gvf_chromosome_counter),
        vcf_chromosome_count=dict(vcf_chromosome_counter) # after merging
    )
    return conversion_payload

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    parser.add_argument("-a", "--assembly", help="FASTA assembly file")
    parser.add_argument("--log", help="Path to log file")
    args = parser.parse_args()

    # Set up logging functionality
    if args.log:
        log_cfg.add_file_handler(args.log)
        logger.info(f"The log file is {args.log}")
    else:
        log_cfg.add_stdout_handler()
    convert(args.gvf_input, args.vcf_output, args.assembly)



if __name__ == "__main__":
    main()
