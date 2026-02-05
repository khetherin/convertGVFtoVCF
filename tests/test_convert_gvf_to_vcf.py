import os.path
import unittest


from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_unstructured_line, \
    convert_gvf_pragmas_for_vcf_header, generate_vcf_header_line, parse_pragma, get_pragma_name_and_value, \
    get_pragma_tokens, \
    get_sample_name_from_pragma, get_unique_sample_names, convert_gvf_pragmas_to_vcf_header, \
    convert_gvf_pragma_comment_to_vcf_header, generate_vcf_header_structured_lines, convert
from convert_gvf_to_vcf.utils import read_in_gvf_header


class TestConvertGVFtoVCF(unittest.TestCase):
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
        self.ordered_list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']

    # the tests below relate to the VCF header (Part 1)
    def test_generate_vcf_header_structured_lines(self):
        all_possible_lines = generate_vcf_header_structured_lines("INFO", self.reference_lookup.mapping_attribute_dict)
        assert isinstance(all_possible_lines, dict)
        assert all_possible_lines["AA"] == '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">'

    def test_generate_vcf_header_unstructured_lines(self):
        key_to_test = "fileformat"
        value_to_test= "VCFv4.4"
        actual_output = generate_vcf_header_unstructured_line(key_to_test, value_to_test)
        assert actual_output == "##fileformat=VCFv4.4"

    def test_get_sample_names_from_pragma(self):
        pragma_name = "#sample"
        pragma_value = "sample_name=JenMale6;subject_name=JenMale6"
        list_of_samples = get_sample_name_from_pragma(pragma_name, pragma_value)
        assert isinstance(list_of_samples, str)
        expected_value = 'JenMale6'
        assert list_of_samples == expected_value

    def test_get_unique_sample_names(self):
        duplicated_samples = ["JenMale6", "JenMale6", "JenMale7"]
        unique_samples = get_unique_sample_names(duplicated_samples)
        expected_list = ["JenMale6", "JenMale7"]
        assert isinstance(unique_samples, list)
        assert unique_samples == expected_list

    def test_convert_gvf_pragmas_to_vcf_header(self):
        list_of_gvf_pragmas_to_convert = [
            '##gff-version 3',
            '##gvf-version 1.06',
            '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955',
            '##file-date 2015-07-15',
            '##genome-build NCBI GRCz10'
        ]
        list_of_gvf_pragmas = ['##file-date', '##gff-version', '##gvf-version', '##species', '##genome-build']
        list_of_converted_pragmas = convert_gvf_pragmas_to_vcf_header(
            list_of_gvf_pragmas_to_convert,
            list_of_gvf_pragmas,
            self.reference_lookup.pragma_to_vcf_map
                                  )
        assert isinstance(list_of_converted_pragmas, list)
        expected_list = ['##gff-version=3', '##gvf-version=1.06', '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10']
        assert list_of_converted_pragmas == expected_list

    def test_convert_gvf_pragma_comment_to_vcf_header(self):
        gvf_pragma_comments_to_convert = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        list_of_gvf_pragma_comments = ['#sample', '#Study_accession', '#Study_type', '#Display_name', '#Publication#Study', '#Assembly_name',
         '#subject']
        list_of_converted_pragma_comments, sample_names_from_pragma_comments = convert_gvf_pragma_comment_to_vcf_header(gvf_pragma_comments_to_convert,
                                                                                           list_of_gvf_pragma_comments,
                                                                                           self.reference_lookup.pragma_to_vcf_map)
        expected_list_of_converted_pragma_comments = ['##Study_accession=nstd62', '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992', '##Journal=Proceedings of the National Academy of Sciences of the United States of America', '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.', '##Publication_year=2012', '##First_author=Kim Brown', '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3', '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']
        expected_list_of_sample_names = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        assert isinstance(list_of_converted_pragma_comments, list)
        assert isinstance(sample_names_from_pragma_comments, list)
        assert list_of_converted_pragma_comments == expected_list_of_converted_pragma_comments
        assert sample_names_from_pragma_comments == expected_list_of_sample_names

    def test_convert_gvf_pragmas_for_vcf_header(self):
        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06', '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        gvf_pragma_comments = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        unique_converted_pragmas, unique_sample_name =convert_gvf_pragmas_for_vcf_header(gvf_pragma, gvf_pragma_comments, self.reference_lookup)
        expected_unique_converted_pragmas = ['##fileformat=VCFv4.4', '##gff-version=3', '##gvf-version=1.06', '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15', '##genome-build=NCBIGRCz10', '##Study_accession=nstd62', '##Study_type=Control Set', '##Display_name=Brown_et_al_2012', '##PMID=22203992', '##Journal=Proceedings of the National Academy of Sciences of the United States of America', '##Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.', '##Publication_year=2012', '##First_author=Kim Brown', '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male', '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3', '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']
        expected_unique_sample_name = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        # test for mandatory first line
        expected_first_line = "##fileformat=VCFv4.4"
        assert unique_converted_pragmas[0] == expected_first_line
        for converted_pragma in unique_converted_pragmas:
            self.assertTrue(converted_pragma in expected_unique_converted_pragmas)

        # assert unique_converted_pragmas == expected_unique_converted_pragmas
        assert unique_sample_name == expected_unique_sample_name

    # the test below relates to the VCF headerline (Part 2)
    def test_generate_vcf_header_line(self):
        header_fields = generate_vcf_header_line(False, self.ordered_list_of_samples)
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'

        header_fields = generate_vcf_header_line(True, self.ordered_list_of_samples)
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'


    # the tests below relate to the GVF header
    def test_parse_pragma(self):
        # testing: pragma has a name and value
        gvf_pragma = "##file-date 2015-07-15"
        delimiter = " "
        name, value = parse_pragma(gvf_pragma, delimiter)
        assert name == "##file-date"
        assert value == "2015-07-15"
        # testing: pragma has only a name, no value, should print warning
        name_only_pragma = "##file-date"
        name, value = parse_pragma(name_only_pragma, delimiter)
        assert name == "##file-date"
        assert value is None
        # testing: invalid pragmas
        invalid_pragma = None
        with self.assertRaises(AttributeError):
            parse_pragma(invalid_pragma, delimiter)

    def test_get_pragma_name_and_value(self):
        pragma_to_test = "##file-date 2015-07-15"
        delimiter = " "
        list_of_pragma = ["##file-date", "##gff-version", "##gvf-version", "##species", "##genome-build"]
        pragma_to_vcf_map = {'##file-date': 'fileDate', '##gff-version': 'gff-version', '##gvf-version': 'gvf-version', '##species': 'species', '##genome-build': 'genome-build', '#sample': 'sample', '#Study_accession': 'Study_accession', '#Study_type': 'Study_type', '#Display_name': 'Display_name', '#Publication': 'Publication', '#Study': 'Study', '#Assembly_name': 'Assembly_name', '#subject': 'subject'}
        vcf_header_key, pragma_name, pragma_value = get_pragma_name_and_value(pragma_to_test, delimiter, list_of_pragma, pragma_to_vcf_map)
        assert vcf_header_key == "fileDate"
        assert pragma_name == "##file-date"
        assert pragma_value == "2015-07-15"

    def test_get_pragma_tokens(self):
        pragma_value = "First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants"
        pragma_tokens = get_pragma_tokens(pragma_value, ";", "=")
        assert len(pragma_tokens) == 2
        # Testing: expected
        assert pragma_tokens[0][0] == "First_author"
        assert pragma_tokens[0][1] == "Kim Brown"
        assert pragma_tokens[1][0] == "Description"
        assert pragma_tokens[1][1] == "Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants"
        # Testing: not expected
        unexpected_pragma_tokens = [['A', '1'], ['B', '2']]
        with self.assertRaises(AssertionError):
            self.assertEqual(unexpected_pragma_tokens, pragma_tokens)

    def test_convert(self):
        # convert_payload = convert(self.input_file , self.output_file, self.assembly)
        # print(convert_payload)
        #####VCF#####
        header_lines = []
        data_lines = []
        with open(self.output_file) as open_file:
            for line in open_file:
                if line.startswith('#'):
                    header_lines.append(line.strip())
                else:
                    data_lines.append(line.strip())
        assert header_lines[0] == '##fileformat=VCFv4.4'
        assert header_lines[-1] == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'
        assert data_lines[0].split('\t')[4] == '<DEL>'
        #####GVF#####
        gvf_header = []
        gvf_features = []
        with open(self.input_file) as input:
            for gvfline in input:
                if gvfline.startswith("#"):
                    gvf_header.append(line.strip())
                else:
                    gvf_features.append(line.strip())
        # checking payload
        # assert convert_payload.samples == ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        # assert convert_payload.vcf_header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'
        # assert convert_payload.gvf_line_count == len(gvf_features)
        # assert convert_payload.vcf_line_count == len(data_lines)


if __name__ == '__main__':
    unittest.main()
