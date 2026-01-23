import os
import unittest


from convert_gvf_to_vcf.utils import read_yaml, read_pragma_mapper, generate_symbolic_allele_dict, \
    build_iupac_ambiguity_code, read_in_gvf_header, read_in_gvf_data
from convert_gvf_to_vcf.lookup import Lookup

class TestUtils(unittest.TestCase):
    def setUp(self):
        # Prepare Directories
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        self.etc_folder =  os.path.join(self.input_folder_parent, "etc")
        input_folder = os.path.dirname(__file__)
        # Prepare Inputs
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # Prepare Outputs
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        # Prepare References
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly)


    def test_read_yaml(self):
        test_yaml_dictionary = read_yaml(os.path.join(self.etc_folder, 'attribute_mapper.yaml'))
        assert len(test_yaml_dictionary) > 0

    def test_read_pragma_mapper(self):
        pragma_to_vcf_header = read_pragma_mapper(os.path.join(self.etc_folder, 'pragma_mapper.tsv'))
        assert len(pragma_to_vcf_header) > 0

    def test_generate_symbolic_allele_dict(self):
        symbolic_allele_dictionary = generate_symbolic_allele_dict(self.reference_lookup.mapping_attribute_dict)
        assert len(symbolic_allele_dictionary) > 0

    def test_read_in_gvf_header(self):
        gvf_pragmas, gvf_non_essentials = read_in_gvf_header(self.input_file)
        assert gvf_pragmas[0] == '##gff-version 3'
        assert all((gvf_pragma.startswith('##') for gvf_pragma in gvf_pragmas))
        assert len(gvf_pragmas) == 5
        assert all((gvf_non_essential.startswith('#') for gvf_non_essential in gvf_non_essentials))
        assert len(gvf_non_essentials) == 15

    def test_read_in_gvf_data(self):
        gvf_features_gen = read_in_gvf_data(self.input_file)
        assert type(gvf_features_gen).__name__ ==  'generator'
        # 7 lines in the GVF
        lines = list(gvf_features_gen)
        assert len(lines) == 7
        assert type(lines[0]).__name__ ==  'GvfFeatureline'
        assert lines[0].start == '1'
        assert lines[0].end == '2'

    def test_build_iupac_ambiguity_code(self):
        expected_dictionary_iupac ={
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'S': ['C', 'G'],
            'D': ['A', 'G', 'T'],
            'W': ['A', 'T'],
            'H': ['A', 'C', 'T'],
            'B': ['C', 'G', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T']
        }
        dictionary_iupac = build_iupac_ambiguity_code()
        assert dictionary_iupac == expected_dictionary_iupac