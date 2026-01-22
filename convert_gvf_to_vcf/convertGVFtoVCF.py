import argparse
import os

from convert_gvf_to_vcf.utils import read_in_gvf_file
from convert_gvf_to_vcf.vcfline import VcfLineBuilder
from convert_gvf_to_vcf.logger import set_up_logging, logger
from convert_gvf_to_vcf.lookup import Lookup
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
            sample_names_from_pragma_comments.append(get_sample_name_from_pragma(pragma_name, pragma_value))
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
    list_of_gvf_pragma_comments = ["#sample", "#Study_accession", "#Study_type", "#Display_name", "#Publication"
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

# This is the main conversion logic
def convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, reference_lookup, ordered_list_of_samples):
    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param reference_lookup: an object that stores important dictionaries to be used for reference lookups.
    :return: standard_header_lines,  list_of_vcf_objects: header lines for this VCF, datalines for this VCF and a list of VCF objects
    """
    list_of_vcf_objects = []
    # Create data structure to store the header lines for this VCF file (standard meta-information lines)
    standard_header_lines ={
        "ALT": [],
        "INFO": [],
        "FILTER": [],
        "FORMAT": [],
    }
    #TODO: place the all_header_lines_per_type_dict into the reference_lookup.

    # Create data structure to store all possible outcomes for header lines (for fields ALT, INFO, FILTER, FORMAT)
    all_header_lines_per_type_dict = {
        htype: generate_vcf_header_structured_lines(htype, reference_lookup.mapping_attribute_dict) for htype in ["ALT", "INFO", "FILTER", "FORMAT"]
    }
    vcf_builder = VcfLineBuilder(standard_header_lines, all_header_lines_per_type_dict, reference_lookup, ordered_list_of_samples)
    # Create a vcf object for every feature line in the GVF (1:1)
    for gvf_featureline in gvf_lines_obj_list:
        vcf_object = vcf_builder.build_vcf_line(gvf_featureline)
        standard_header_lines = vcf_builder.build_vcf_header() # Forms headerlines from vcf line obj
        # Store VCF object in the list
        list_of_vcf_objects.append(vcf_object)
    # Returns the header of the VCF file, and the objects.
    return standard_header_lines, list_of_vcf_objects

# The functions below relate to the VCF objects
def collect_missing_format_flags(list_of_vcf_objects):
    """ Returns a list of booleans for each VCF object (True if format is missing, False if format keys present)
    :params: list_of_vcf_objects: list of vcf objects
    :return: missing_format_flags: list of booleans (True if format is missing, False if format keys present)
    """
    missing_format_flags = []
    for index in range(1, len(list_of_vcf_objects)):
        if list_of_vcf_objects[index].format_keys == ['.']:
            format_flag = True
            missing_format_flags.append(format_flag)
        else:
            format_flag = False
            missing_format_flags.append(format_flag)
    return missing_format_flags


def compare_vcf_objects(list_of_vcf_objects):
    """ Compares VCF objects in the list with the VCF object before it. Returns boolean values.
    :params: list_of_vcf_objects: list of vcf objects
    :return: comparison_results: list of booleans. For future reference, if True, this will determine merging lines; if False, this will determine use of the previous line.
    """
    comparison_results = []
    # For each vcf line object, compare with the previous vcf line object in the list
    for index in range(1, len(list_of_vcf_objects)):
        current_vcf_object = list_of_vcf_objects[index]
        previous_vcf_object = list_of_vcf_objects[index - 1]
        # Determines the VCF line objects as equal based on the CHROM, POS and REF being the same (__eq__ in Vcfline)
        if current_vcf_object == previous_vcf_object:
            comparison_results.append(True) # This will use require merging.
        else:
            comparison_results.append(False) # No merging required. Use previous object.
    return comparison_results

def merge_vcf_objects(previous, current, list_of_sample_names):
    """ Merge VCF objects.
    :params: previous: previous VCF line object
    :params: current: current VCF line object
    :params: list_of_sample_names: sample names
    :return: merged_object
    """
    merged_object = previous.merge(current, list_of_sample_names)
    return merged_object

def keep_vcf_objects(previous, list_of_sample_names):
    """ Keep VCF objects.
    :params: previous VCF line object
    :return: kept_object
    """
    kept_object = previous.keep(list_of_sample_names)
    return kept_object

def determine_merge_or_keep_vcf_objects(list_of_vcf_objects, comparison_results, list_of_sample_names):
    """ Runs through the list of VCF objects and its corresponding comparison result.
    If True, merge parts of the vcf object together. If False, use the previous object
    :params: list_of_vcf_objects: list of vcf line objects
    :return: merge_or_kept_objects: list of vcf line objects that have either been merged or kept as is.
    """
    merge_or_kept_objects = []
    # start at 1 to ensure the first element has a previous object
    for index, compare_result in enumerate(comparison_results, start=1):
        # Merge if the previous and current VCF object are the same (compare_result is True)
        if compare_result:
            merged_object = merge_vcf_objects(list_of_vcf_objects[index - 1], list_of_vcf_objects[index], list_of_sample_names)
        # Keep previous if previous and current VCF object are different (compare_result is False)
        else:
                # keep the previous VCF line object
            kept_object = keep_vcf_objects(list_of_vcf_objects[index - 1], list_of_sample_names)
            merge_or_kept_objects.append(kept_object)
    merge_or_kept_objects.append(list_of_vcf_objects[-1])
    return merge_or_kept_objects

def get_chrom_pos_of_vcf_object(obj):
    """
    Returns chromosome and position of an object. Used to help sort VCF object by chromosome name and by numeric position.
    :params: obj : vcf line object
    :return: obj.chrom, obj.pos
    """
    return (obj.chrom, int(obj.pos))

def has_duplicates(list_of_objects):
    """ Checks chrom and pos of vcf line objects to see if there are duplicates. If list of vcf object has duplicates, merge again.
    :params: list_of_objects
    :return: duplicate_flag - boolean value (True = has duplicates so merge again, False= no duplicates)
    """
    duplicate_flag = False
    list_of_duplicate_chrom_pos = []
    seen_chrom_pos = set()
    for obj in list_of_objects:
        chrom_pos = (obj.chrom, int(obj.pos))
        if chrom_pos in seen_chrom_pos:
            # returns true to allow to merge again
            duplicate_flag = True
            list_of_duplicate_chrom_pos.append(chrom_pos)
        else:
            seen_chrom_pos.add(chrom_pos)
    return duplicate_flag, list_of_duplicate_chrom_pos


def get_list_of_merged_vcf_objects(list_of_vcf_objects, samples):
    """ Compares VCF objects, merges VCF objects, sorts VCF objects. This gives a list of sorted merged vcf objects.
    :params: list_of_vcf_objects
    :params: samples
    :returns: merge_or_kept_vcf_objects
    """
    comparison_flags = compare_vcf_objects(list_of_vcf_objects)  # Identifies which VCF objects to merge
    merge_or_kept_vcf_objects = determine_merge_or_keep_vcf_objects(list_of_vcf_objects, comparison_flags, samples)
    merge_or_kept_vcf_objects.sort(key=get_chrom_pos_of_vcf_object)  # sorting by chromosome and position
    return merge_or_kept_vcf_objects

def filter_duplicates_by_merging(chrom_pos_list, has_dups, list_of_vcf_objects,
                                 list_of_vcf_objects_to_be_filtered, samples):
    """
    :params: chrom_pos_list: list of tuples - this represents duplicate positions
    :params: has_dups: boolean - True if it contains duplicates
    :params: list_of_vcf_objects: list of VCF objects (GVF converted to VCF; with no merging/remove dups)
    :params: list_of_vcf_objects_to_be_filtered: list of VCF objects from a previous merge
    :params: samples: names
    :returns: filtered_merge_or_kept_vcf_objects - list of vcf objects with duplicates from this iteration removed
    """
    if has_dups:
        for chrom_pos in chrom_pos_list:
            chrom_to_search = chrom_pos[0]
            pos_to_search = chrom_pos[1]
            vcf_objects_to_merge = []
            for vcf_object in list_of_vcf_objects:
                if vcf_object.chrom == chrom_to_search and vcf_object.pos == pos_to_search:
                    vcf_objects_to_merge.append(vcf_object)

    merge_duplicates = get_list_of_merged_vcf_objects(vcf_objects_to_merge, samples)
    filtered_merge_or_kept_vcf_objects = [x for x in list_of_vcf_objects_to_be_filtered if x not in vcf_objects_to_merge]
    filtered_merge_or_kept_vcf_objects.extend(merge_duplicates)
    filtered_merge_or_kept_vcf_objects.sort(key=get_chrom_pos_of_vcf_object)
    return filtered_merge_or_kept_vcf_objects

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
        log_path = set_up_logging(args.log)
    else:
        log_path = set_up_logging()

    # Log the inputs and outputs.
    logger.info("Running the GVF to VCF converter")
    logger.info(f"The provided input file is: {args.gvf_input}")
    logger.info(f"The provided output file is: {args.vcf_output}")
    if args.assembly:
        logger.info(f"The provided assembly file is: {args.assembly}")
    assembly_file = os.path.abspath(args.assembly)
    assert os.path.isfile(assembly_file), "Assembly file does not exist"
    logger.info(f"The log file is {log_path}")

    # Read input file and separate out its components
    logger.info(f"Reading in the following GVF input: {args.gvf_input}")
    gvf_pragmas, gvf_pragma_comments, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)
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

    # Convert each feature line in the GVF file to a VCF object (stores all the data for a line in the VCF file).
    # NOTE: Main Logic lives here.
    (header_lines_per_type,list_of_vcf_objects) = convert_gvf_features_to_vcf_objects(gvf_lines_obj_list, reference_lookup, ordered_list_of_samples=samples)
    logger.info(f"Writing to the following VCF output: {args.vcf_output}")
    logger.info("Generating the VCF header and the meta-information lines")
    with open(args.vcf_output, "w") as vcf_output:

        logger.info(f"Total number of samples in this VCF: {len(samples)}")

        # Part 1 of VCF file: Write the VCF header. This will include preserved data from the GVF file.
        for pragma in pragmas_for_vcf:
            vcf_output.write(f"{pragma}\n")
        for header_type in header_lines_per_type.values():
            for header_line in header_type:
                vcf_output.write(f"{header_line}\n")

        # Part 2 of VCF file: Write the VCF header line.
        # Determine if the header is the 8 mandatory fields or 8 mandatory fields + FORMAT + sample names.
        missing_flags = collect_missing_format_flags(list_of_vcf_objects) # True if format keys are missing, False if present
        if all(missing_flags):
            is_missing_format_value = True
            logger.info("No Format Keys detected. Printing mandatory VCF headers.")
        else:
            is_missing_format_value = False
        # Write the header.
        header_fields = generate_vcf_header_line(is_missing_format=is_missing_format_value, samples=samples)
        vcf_output.write(f"{header_fields}\n")

        # Part 3 of VCF file: Write the VCF data lines. This will contain info about the position in the genome,
        # its variants and genotype information per sample.
        if (list_of_vcf_objects):
            logger.info("Generating the VCF datalines")
            # Each GVF feature has been converted to a VCF object so begin comparing and merging the VCF objects.
            # initial merge
            merge_or_kept_vcf_objects = get_list_of_merged_vcf_objects(list_of_vcf_objects, samples)
            # identify if duplicates are present after merging
            has_dups, chrom_pos_list = has_duplicates(merge_or_kept_vcf_objects)
            # while duplicates are present, merge, then re-check for dups
            max_iterations = 100
            iteration = 0
            list_of_vcf_objects_to_be_filtered = merge_or_kept_vcf_objects
            while has_dups and iteration < max_iterations:
                filtered_merge_or_kept_vcf_objects = filter_duplicates_by_merging(chrom_pos_list, has_dups,
                                                                                  list_of_vcf_objects,
                                                                                  list_of_vcf_objects_to_be_filtered, samples)
                has_dups, chrom_pos_list = has_duplicates(filtered_merge_or_kept_vcf_objects)
                iteration += 1
                list_of_vcf_objects_to_be_filtered = filtered_merge_or_kept_vcf_objects
                logger.info(f"Iteration of merge (remove dups): {iteration}")

            # Write the VCF objects as data lines in the VCF file.
            if iteration != 0:
                for vcf_line_object in filtered_merge_or_kept_vcf_objects:
                    vcf_output.write(str(vcf_line_object) + "\n")
            else:
                for vcf_line_object in merge_or_kept_vcf_objects:
                    vcf_output.write(str(vcf_line_object) + "\n")
        else:
            logger.warning("No feature lines were found for this GVF file.")
    vcf_output.close()
    logger.info("GVF to VCF conversion complete")

if __name__ == "__main__":
    main()
