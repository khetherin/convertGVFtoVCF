import argparse
import os

from Bio import SeqIO
from convert_gvf_to_vcf.utils import read_file, read_info_attributes, read_pragma_mapper, read_sequence_ontology_symbolic_allele, read_in_gvf_file
from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.assistingconverter import generate_custom_structured_meta_line
# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')

def generate_vcf_header_structured_lines(header_type):
    """ Generates a dictionary of all possible standard structured lines for INFO/FILTER/FORMAT/ALT
    :param header_type: type of header file to read i.e. ALT, FILTER, INFO or FORMAT
    :return: dictionary of all possible standard structured lines keys for the header type
    """
    all_possible_lines = {}
    prefix = {}
    prefix["reserved"] = True
    prefix["sv"] = True
    if header_type == 'ALT':
        prefix["reserved"] = False
    if prefix["reserved"]:
        reserved_key = read_file("reserved", header_type)
        for r_key in reserved_key:
            key_id = reserved_key[r_key][0]
            number = reserved_key[r_key][1]
            type_for_key = reserved_key[r_key][2]
            description = reserved_key[r_key][3]
            reserved_string = f'##{header_type}=<ID={key_id},Number={number},Type={type_for_key},Description="{description}">'
            all_possible_lines[key_id] = reserved_string
    if prefix['sv']:
        sv_key = read_file("sv", header_type)
        for s_key in sv_key:
            sv_key_id = sv_key[s_key][0]
            sv_line = sv_key[s_key][1]
            all_possible_lines[sv_key_id] = sv_line
    return all_possible_lines

def generate_custom_unstructured_meta_line(vcf_unstructured_key,
                                           vcf_unstructured_value):
    """ Generates a formatted unstructured metainformation line using a custom key value pair.
    :param vcf_unstructured_key: key for custom unstructured metainformation line
    :param vcf_unstructured_value: value for custom unstructured metainformation line
    :return: custom_unstructured_string
    """
    custom_unstructured_string = f"##{vcf_unstructured_key}={vcf_unstructured_value}"
    return custom_unstructured_string


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
            print("WARNING: no value for the following pragma", pragma_to_parse)
        return pragma_name, pragma_value
    except ValueError:
        print("Skipping this, can't be parsed", pragma_to_parse)

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
        pragma_tokens = element.split(second_delimiter)
    return pragma_tokens

#step 9 using custom unstructured meta-information line = generate_custom_unstructured_metainfomation_line
def generate_vcf_metainformation(gvf_pragmas, gvf_non_essential, list_of_vcf_objects,
                                 standard_lines_dictionary):
    """ Generates a list of metainformation lines for the VCF header
    :param gvf_pragmas: list of gvf pragmas to convert
    :param gvf_non_essential: list of non-essential gvf pragmas to convert
    :param list_of_vcf_objects: list of vcf objects
    :param standard_lines_dictionary: dictionary of standard lines
    :return: unique_pragmas_to_add, sample_names: a list of pragmas (this list contains no duplicates), list of sample names
    """
    pragmas_to_add = []
    unique_pragmas_to_add = []
    sample_names = []
    unique_alt_lines_to_add = []
    unique_info_lines_to_add = []
    unique_filter_lines_to_add = []
    unique_format_lines_to_add = []
    # MANDATORY: file format for VCF
    pragmas_to_add.append(generate_custom_unstructured_meta_line("fileformat", "VCFv4.4"))
    #Go through essential pragmas
    #TODO: list of pragmas to add:reference=file, contig, phasing,INFO#
    list_of_pragma = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
    pragma_to_vcf_map = read_pragma_mapper(os.path.join(etc_folder, 'pragma_mapper.tsv'))
    for pragma in gvf_pragmas:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma, " ", list_of_pragma, pragma_to_vcf_map)
        pragmas_to_add.append(generate_custom_unstructured_meta_line(vcf_header_key, pragma_value))
        # adding in source
        for vcf_obj in list_of_vcf_objects:
            pragmas_to_add.append(generate_custom_unstructured_meta_line("source", vcf_obj.source))
    # Go through non-essential pragmas
    list_of_non_essential_pragma = ["#sample", "#Study_accession", "#Study_type", "#Display_name", "#Publication"
                                    "#Study", "#Assembly_name", "#subject"]
    for non_essential_pragma in gvf_non_essential:
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(non_essential_pragma, ": ", list_of_non_essential_pragma, pragma_to_vcf_map)
        if pragma_name.startswith("#Publication"):
            publication_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_custom_unstructured_meta_line(publication_tokens[0], publication_tokens[1]))
        elif pragma_name == "#Study":
            study_tokens = get_pragma_tokens(pragma_value, ";", "=")
            pragmas_to_add.append(generate_custom_unstructured_meta_line(study_tokens[0], study_tokens[1]))
        else:
            if vcf_header_key is not None:
                pragmas_to_add.append(generate_custom_unstructured_meta_line(vcf_header_key, pragma_value))
        # populating sample headers
        if pragma_name.startswith("#sample"):
            list_of_sample_information = pragma_value.split(";")
            for sample_info in list_of_sample_information:
                if sample_info.startswith("sample_name"):
                    sample_names.append(sample_info.split("=")[1])

    print("Total number of samples in this VCF: ", len(sample_names))

    unique_pragmas_to_add = list(dict.fromkeys(pragma for pragma in pragmas_to_add if pragma not in unique_pragmas_to_add))
    unique_alt_lines_to_add = list(dict.fromkeys(alt_line for alt_line in standard_lines_dictionary["ALT"] if alt_line not in unique_alt_lines_to_add))
    unique_info_lines_to_add = list(dict.fromkeys(info_line for info_line in standard_lines_dictionary["INFO"] if info_line not in unique_info_lines_to_add))
    unique_filter_lines_to_add = list(dict.fromkeys(filter_line for filter_line in standard_lines_dictionary["FILTER"] if filter_line not in unique_filter_lines_to_add))
    unique_format_lines_to_add = list(dict.fromkeys(format_line for format_line in standard_lines_dictionary["FORMAT"] if format_line not in unique_format_lines_to_add))

    return unique_pragmas_to_add, sample_names, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add

# step 10
def generate_vcf_header_line(samples):
    """ Generates the VCF header line
    :param samples: list of samples, these will appear in the header line
    :return: vcf_header: a string
    """
    vcf_header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for sample in samples:
        vcf_header_fields.append(sample)
    vcf_header = '\t'.join(vcf_header_fields)
    return vcf_header

def gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                assembly_file):

    """ Creates VCF objects from GVF feature lines and stores the VCF objects.
    :param gvf_lines_obj_list: list of GVF feature line objects
    :param assembly_file: FASTA file to assembly
    :return: header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects: dictionary of lists and a list of VCF objects
    """
    vcf_data_lines = {}  # DICTIONARY OF LISTS
    list_of_vcf_objects = []
    # standard meta-information lines for this VCF file
    header_standard_lines_dictionary ={
        "ALT": [],
        "INFO": [],
        "FILTER": [],
        "FORMAT": [],
    }
    all_header_lines_per_type_dict = {
        htype: generate_vcf_header_structured_lines(htype) for htype in ["ALT", "INFO", "FILTER", "FORMAT"]
    }
    info_attribute_dict = read_info_attributes(os.path.join(etc_folder, 'INFOattributes.tsv'))  # needed to generate custom strings

    symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(os.path.join(etc_folder, 'svALTkeys.tsv'))

    # create a vcf object for every feature line in the GVF (1:1)
    # add the newly created vcf object to the vcf data line it belongs to
    # (1:many; key=chrom_pos; 1 key: many vcf objects)
    for gvf_featureline in gvf_lines_obj_list:
        vcf_object = VcfLine(gvf_featureline,
                             info_attribute_dict,
                             symbolic_allele_dictionary,
                             assembly_file,
                             header_standard_lines_dictionary,
                             all_header_lines_per_type_dict)

        list_of_vcf_objects.append(vcf_object)
        if vcf_object.key in vcf_data_lines:
            vcf_data_lines[vcf_object.key].append(vcf_object)
        else:
            vcf_data_line_objects_list = []
            vcf_data_line_objects_list.append(vcf_object)
            vcf_data_lines[vcf_object.key] = vcf_data_line_objects_list
        # check the number of objects to see if they are merged
        # for key in vcf_data_lines.keys():
        #     vcf_obj_list = vcf_data_lines[key]
            # print("for ", key, " the number of vcf objects is: ", len(vcf_obj_list))
    return header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects


def populate_sample_formats(list_of_sample_names):
    """ Populates a dictionary using a list of sample names. Dictionary key is sample name, value is the sample's format value.
    :param list_of_sample_names: list of sample names
    :return:sample_name_format_value: dictionary of sample names => sample format value
    """
    sample_name_format_value = {}
    for sample in list_of_sample_names:
        sample_name_format_value[sample] = "sampleFORMAThere" #TODO: fill this in
    return sample_name_format_value

def format_sample_values(sample_name_format_value):
    """ Creates a partial vcf data line of sample format values.
    :param sample_name_format_value: dictionary of sample names => sample format value
    :return: sample_format_values_string: formatted string
    """
    sample_format_values_string = ""
    for key in sample_name_format_value:
        sample_format_values_string = sample_format_values_string + sample_name_format_value[key] + "\t"
    return sample_format_values_string

def format_vcf_datalines(list_of_vcf_objects, list_of_sample_names):
    """ Iterates through a list of VCF objects and sample names and formats them as a VCF dataline.
    :param list_of_vcf_objects: list of vcf objects
    :param list_of_sample_names: list of sample names
    :return: formatted_vcf_datalines: list of formatted vcf datalines
    """
    sample_name_format_value = populate_sample_formats(list_of_sample_names)
    sample_format_values_string = format_sample_values(sample_name_format_value)

    formatted_vcf_datalines = []
    for vcf_obj in list_of_vcf_objects:
        vcf_info_string = ";".join(vcf_obj.info)
        vcf_line = (f"{vcf_obj.chrom}\t"
                        f"{vcf_obj.pos}\t"
                        f"{vcf_obj.id}\t"
                        f"{vcf_obj.ref}\t" #TODO: should this always be empty
                        f"{vcf_obj.alt}\t" #TODO: should this always be empty
                        f"{vcf_obj.qual}\t" #TODO: should this always be empty
                        f"{vcf_obj.filter}\t" #TODO: should this always be empty
                        #f"{vcf_obj.info}\t"
                        f"{vcf_info_string}\t"
                        f"{vcf_obj.format}\t"
                        f"{sample_format_values_string}"
                        )
        formatted_vcf_datalines.append(vcf_line)
    return formatted_vcf_datalines

def main():
    print("Running the GVF to VCF converter")
    # step 1
    parser = argparse.ArgumentParser()
    parser.add_argument("gvf_input", help="GVF input file.")
    parser.add_argument("vcf_output", help="VCF output file.")
    parser.add_argument("-a", "--assembly", help="FASTA assembly file")
    args = parser.parse_args()
    print("The provided input file is: ", args.gvf_input)
    print("The provided output file is: ", args.vcf_output)
    if args.assembly:
        print("The provided assembly file is: ", args.assembly)
    assembly_file = os.path.abspath(args.assembly)
    assert os.path.isfile(assembly_file), "Assembly file does not exist"

    # custom meta-information lines for this VCF file
    gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(args.gvf_input)

    (
        header_standard_lines_dictionary,
        vcf_data_lines,
        list_of_vcf_objects
    ) = gvf_features_to_vcf_objects(gvf_lines_obj_list, assembly_file)


    print("Writing to the following VCF output: ", args.vcf_output)
    print("Generating the VCF header and the meta-information lines")

    with open(args.vcf_output, "w") as vcf_output:
        (
            unique_pragmas_to_add,
            samples,
            unique_alt_lines_to_add,
            unique_info_lines_to_add,
            unique_filter_lines_to_add,
            unique_format_lines_to_add
        ) = generate_vcf_metainformation(gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary)
        for pragma in unique_pragmas_to_add:
            vcf_output.write(f"{pragma}\n")
        for alt_lines in unique_alt_lines_to_add:
            vcf_output.write(f"{alt_lines}\n")
        for info_lines in unique_info_lines_to_add:
            vcf_output.write(f"{info_lines}\n")
        for filter_lines in unique_filter_lines_to_add:
            vcf_output.write(f"{filter_lines}\n")
        for format_lines in unique_format_lines_to_add:
            vcf_output.write(f"{format_lines}\n")
        header_fields = generate_vcf_header_line(samples)
        vcf_output.write(f"{header_fields}\n")
        print("Generating the VCF datalines")
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        for line in formatted_vcf_datalines:
            vcf_output.write(f"{line}\n")
    vcf_output.close()
    print("GVF to VCF conversion complete")


if __name__ == "__main__":
    main()
