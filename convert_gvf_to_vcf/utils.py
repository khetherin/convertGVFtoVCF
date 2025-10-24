# this file contains readers only
import os
from convert_gvf_to_vcf.gvffeature import GvfFeatureline
# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')


def read_file(prefix, header_type):
    """Reads in {reserved/sv}{INFO/FORMAT}keys.tsv files and returns the dictionary where the key is the KEYID
    (usually column 1) and the value is a list of file tokens
    :param prefix: prefix of files to read i.e. sv or reserved
    :param header_type: type of header file to read i.e. INFO or FORMAT
    :return: file lines: a dictionary where the key is the KEYID and the value is a list of file tokens
    """
    file_lines = {}
    keys_tsv_file = os.path.join(etc_folder, f'{prefix}{header_type}keys.tsv')
    print(keys_tsv_file)
    try:
        with open(keys_tsv_file, encoding='utf-8', errors='replace') as keys_file:
            next(keys_file)  # Skip the header
            for line in keys_file:
                file_tokens = line.rstrip().split("\t")
                key_id = file_tokens[0]
                number_of_tokens = len(file_tokens)
                if number_of_tokens >= 2:
                    for token in range(number_of_tokens):
                        value_to_add = file_tokens[token]
                        file_lines.setdefault(key_id, []).append(value_to_add)
    except FileNotFoundError as e:
        print(f'File not found: {keys_tsv_file}')
        raise e
    return file_lines

def read_info_attributes(info_attributes_file):
    """ Read in the file containing specific INFO attributes.
    :param info_attributes_file: A file containing the specific attributes
    :return: attribute_dict: A dictionary of key id and the list of attribute tokens
    """
    attribute_dict = {}  # dictionary of dgva specific INFO attributes
    with open(info_attributes_file) as open_file:
        next(open_file)
        for line in open_file:
            attribute_tokens = line.rstrip().split("\t")
            key = attribute_tokens[0]
            attribute_dict[key] = attribute_tokens
    return attribute_dict

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

def read_sequence_ontology_symbolic_allele(so_symbolic_allele_file):
    """ Read in the file containing sequence ontology symbolic allele and returns a dictionary.
    :param: so_symbolic_allele_file - the file of sequence ontology symbolic alleles.
    :return: symbolic allele dictionary - symbolic alleles as key and list of variant types as the value.
    """
    symbolic_allele_dict = {}
    with open(so_symbolic_allele_file) as so_symbolic_allele:
        next(so_symbolic_allele)
        for line in so_symbolic_allele:
            allele_tokens = line.rstrip().split("\t")
            symb_allele = allele_tokens[0]
            sequence_ontology_id = allele_tokens[2]
            name = allele_tokens[3]
            description = allele_tokens[4]
            symbolic_allele_dict.setdefault(name, []).append(sequence_ontology_id)
            symbolic_allele_dict.setdefault(name, []).append(symb_allele)
            symbolic_allele_dict.setdefault(name, []).append(description)
    return symbolic_allele_dict

def read_in_gvf_file(gvf_input):
    """ Reads in the user provided GVF file.
    :param gvf_input: arguments.gvf_input
    :return: gvf_pragmas, gvf_non_essential, gvf_lines_obj_list
    """
    gvf_pragmas = []  # list of pragma lines starting with: ##
    gvf_non_essential = []  # list of non-essential lines starting with: #
    features = []
    gvf_lines_obj_list = []  # list of objects when reading in gvf files, one object represents a gvf line
    print("Reading in the following GVF input: " + gvf_input)
    with open(gvf_input) as gvf_file:
        for line in gvf_file:
            if line.startswith("##"):
                gvf_pragmas.append(line.rstrip())
            elif line.startswith("#"):
                gvf_non_essential.append(line.rstrip())
            else:
                features.append(line.rstrip())
    for feature in features:
        f_list = feature.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_lines_obj_list.append(line_object)
    return gvf_pragmas, gvf_non_essential, gvf_lines_obj_list

