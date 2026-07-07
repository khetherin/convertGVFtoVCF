import os.path
import unittest


from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.convert_gvf_to_vcf_logic import generate_vcf_header_unstructured_line, \
    convert_gvf_pragmas_for_vcf_header, generate_vcf_header_line, parse_pragma, get_pragma_name_and_value, \
    get_pragma_tokens, \
    get_sample_name_from_pragma, get_unique_sample_names, convert_gvf_pragmas_to_vcf_header, \
    convert_gvf_pragma_comment_to_vcf_header, generate_vcf_header_structured_lines, convert
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestConvertGVFtoVCF(unittest.TestCase):
    def setUp(self):
        # Prepare Directories
        self.paths = ProjectPaths()
        self.etc_folder = self.paths.etc_dir
        self.input_folder = self.paths.package_dir
        self.tests_folder = self.paths.test_dir
        # Prepare Inputs
        self.input_file = os.path.join(self.tests_folder, "input", "zebrafish.gvf")
        # Prepare Outputs
        self.output_file = os.path.join(self.tests_folder, "output", "a.vcf")
        # Prepare References
        self.assembly = os.path.join(self.tests_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly, self.paths)
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
        convert(self.input_file , self.output_file, self.assembly, self.paths)
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
                    gvf_header.append(gvfline.strip())
                else:
                    gvf_features.append(gvfline.strip())
        assert gvf_header[1] == "##gvf-version 1.06"
        assert gvf_features[0] == "chromosome1	DGVa	copy_number_loss	1	2	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        # check statistics file exists and what is inside it
        vcf_output_directory = os.path.dirname(os.path.abspath(self.output_file))
        gvf_filename_only = os.path.basename(self.input_file)
        gvf_file_base, _ = os.path.splitext(gvf_filename_only)
        stats_file = os.path.join(vcf_output_directory, f"{gvf_file_base}.gvf.stats.txt")


        assert os.path.exists(stats_file) == True
        stats_lines = []
        with open(stats_file, 'r') as file:
            for stats_line in file:
                stats_lines.append(stats_line.strip())
        assert stats_lines[12] == "Counter({'chromosome1': 7})"
        assert stats_lines[13] == "Number of GVF feature lines:  	7"
        assert stats_lines[16] == "Counter({'copy_number_gain': 4, 'copy_number_loss': 3})"
        assert stats_lines[20] == "Counter({'chromosome1': 4})"
        assert stats_lines[22] == "Counter({'<DEL>': 3, '<DUP>': 1})"



    def test_convert_500(self):
        input_folder = os.path.join(self.tests_folder, "input")
        output_folder = os.path.join(self.tests_folder, "output")
        # Prepare Inputs
        input_file = os.path.join(input_folder, "drosophila_estd205_lines_500_sorted.gvf")
        output_file = os.path.join(output_folder, "drosophila_estd205_lines_500.vcf")
        os.makedirs(output_folder, exist_ok=True)
        # Prepare References
        assembly = os.path.join(input_folder, "drosophila_GCA_000001215.2_chr4.fa")
        assert os.path.exists(input_file) == True
        assert os.path.exists(assembly) == True
        convert(input_file, output_file, assembly, self.paths)

    def test_convert_position_shift(self):
        """
        Two GVF records at the same GVF position (5) that produce different VCF positions:
          - Record 1 (first in GVF): SNP, Variant_seq=G → VCF pos=5 (no padding)
          - Record 2 (second in GVF): Del, Variant_seq=. → VCF pos=4 (padded before)
        """
        input_folder = os.path.join(self.tests_folder, "input")
        output_folder = os.path.join(self.tests_folder, "output")
        os.makedirs(output_folder, exist_ok=True)
        input_file = os.path.join(input_folder, "zebrafish_position_shift.gvf")
        output_file = os.path.join(output_folder, "zebrafish_position_shift.vcf")

        convert(input_file, output_file, self.assembly, self.paths)

        data_lines = []
        with open(output_file) as f:
            for line in f:
                if not line.startswith('#'):
                    data_lines.append(line.strip().split('\t'))
        assert len(data_lines) == 2, f"Expected 2 VCF data lines, got {len(data_lines)}"
        # After sorting: pos=4 (DEL) must come before pos=5 (gain)
        assert data_lines[0][1] == '4'
        assert data_lines[0][3] == 'TACGT'
        assert data_lines[0][4] == '<DEL>'
        assert data_lines[1][1] == '5'
        assert data_lines[1][3] == 'A'
        assert data_lines[1][4] == 'G'

if __name__ == '__main__':
    unittest.main()
