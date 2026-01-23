import os

from Bio import SeqIO

from convert_gvf_to_vcf.utils import read_yaml, generate_symbolic_allele_dict, build_iupac_ambiguity_code, \
    read_pragma_mapper

# setting up paths to useful directories
convert_gvf_to_vcf_folder = os.path.dirname(__file__)
etc_folder = os.path.join(convert_gvf_to_vcf_folder, 'etc')


class Lookup:
    """
    The class is responsible for the storage of look up dictionaries or files in the etc folder for a VCF file.
    """
    def __init__(self, assembly_file):
        self.mapping_attribute_dict = read_yaml(os.path.join(etc_folder, "attribute_mapper.yaml"))
        self.symbolic_allele_dictionary = generate_symbolic_allele_dict(self.mapping_attribute_dict)
        self.assembly_file = assembly_file
        # using .index for memory efficiency:  https://biopython.org/docs/1.76/api/Bio.SeqIO.html#input-multiple-records
        self.assembly_fasta_indexed = SeqIO.index(self.assembly_file, "fasta")
        self.iupac_ambiguity_dictionary = build_iupac_ambiguity_code()
        self.pragma_to_vcf_map = read_pragma_mapper(os.path.join(etc_folder, 'pragma_mapper.tsv'))
        # self.all_possible_vcf_header_lines_dictionary={
        #         htype: generate_vcf_header_structured_lines(htype, self.mapping_attribute_dict) for htype in ["ALT", "INFO", "FILTER", "FORMAT"]
        #     }
