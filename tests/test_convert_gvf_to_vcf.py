import os.path
import unittest

#from convert_gvf_to_vcf.utils import read_file
from convert_gvf_to_vcf.convertGVFtoVCF import read_file, generate_custom_unstructured_meta_line, read_in_gvf_file, \
    gvf_features_to_vcf_objects, format_vcf_datalines, \
    generate_vcf_metainformation, generate_vcf_header_structured_lines,  \
    generate_vcf_header_line, populate_sample_formats, \
    format_sample_values, read_info_attributes, read_sequence_ontology_symbolic_allele

from convert_gvf_to_vcf.vcfline import VcfLine
from convert_gvf_to_vcf.gvffeature import GvfFeatureline

class TestConvertGVFtoVCF(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'convert_gvf_to_vcf'))
        # the inputs below are INFO attribute files
        self.info_attribute_input_file = os.path.join(self.input_folder_parent, "etc", "INFOattributes.tsv")

        self.symbolic_allele_file = os.path.join(self.input_folder_parent,"etc", 'svALTkeys.tsv')
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")

    def test_read_file(self):
        prefix = "reserved"
        header_type = "FORMAT"
        file_lines = read_file(prefix,header_type)
        assert len(file_lines) > 1

    def test_read_in_gvf_file(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assert len(gvf_pragmas) > 1
        assert len(gvf_non_essential) > 1
        assert len(gvf_lines_obj_list) > 1

    def test_read_info_attributes(self):
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        assert len(info_attribute_dict) > 1
        kv = next(iter(info_attribute_dict.items()))
        key_to_check = next(iter(info_attribute_dict))
        value_to_check = ["3'_inner_flank_link", '.', 'String', 'Three prime inner flank link', 'DGVa']
        assert kv[0] == key_to_check
        assert kv[1] == value_to_check

    def test_read_sequence_ontology_symbolic_allele(self):
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assert len(symbolic_allele_dictionary) > 1

    def test_gvf_features_to_vcf_objects(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        assembly_file = self.assembly

        # standard structured meta-information lines for this VCF file
        header_lines_for_this_vcf, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          assembly_file)
        assert len(vcf_data_lines) > 1
        assert len(list_of_vcf_objects) > 1

    def test_add_padded_base(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly

        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }

        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        test_ref = "A"
        test_alt = "T"
        padded_base, pos, ref, alt = v.add_padded_base(test_ref, test_alt, True)
        assert padded_base is not None
        assert pos is not None
        assert ref is not None
        assert alt is not None

    def test_build_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)

        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        assert len(my_ipuac_dictionary) > 0

    def test_convert_iupac_ambiguity_code(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)


        my_ipuac_dictionary = v.build_iupac_ambiguity_code()
        ref_to_convert = "TAGD"
        converted_ref_allele = v.convert_iupac_ambiguity_code(my_ipuac_dictionary, ref_to_convert)
        assert converted_ref_allele not in ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]

    def test_check_ref(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                     f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)

        reference_allele_to_check = "TGCR"
        new_ref = v.check_ref(reference_allele_to_check)
        iupac_code = ["R", "Y", "M", "K", "S", "D", "W", "H", "B", "V", "D", "N"]
        assert all(code not in new_ref for code in iupac_code)

    def test_get_ref(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        lines_custom_unstructured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        reference_allele = v.get_ref()
        assert len(reference_allele) != 0
        assert reference_allele == 'TA'

    def test_generate_symbolic_allele(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)
        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly
        # custom meta-information lines for this VCF file
        lines_custom_structured = []
        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        output_symbolic_allele, info_field, output_lines_standard_ALT, output_lines_standard_INFO = v.generate_symbolic_allele(standard_lines_dictionary, all_possible_lines_dictionary)
        assert output_symbolic_allele == '<DEL>'
        assert info_field == ['END=81', 'SVLEN=4', 'IMPRECISE', 'CIPOS=0,1', 'CIEND=0,1', 'END=80', 'SVLEN=4', 'IMPRECISE', 'CIPOS=1,2', 'CIEND=1,2']
        assert output_lines_standard_ALT == ['"##ALT=<ID=DEL,Description=""Deletion"">"', '"##ALT=<ID=DEL,Description=""Deletion"">"']
        print("here is ", output_lines_standard_INFO)
        assert output_lines_standard_INFO ==  [
            '##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier.">',
            '##INFO=<ID=Name,Number=.,Type=String,Description="name">',
            '##INFO=<ID=Alias,Number=.,Type=String,Description="A secondary name.">',
            '##INFO=<ID=variant_call_so_id,Number=.,Type=String,Description="variant call SO id">',
            '##INFO=<ID=parent,Number=.,Type=String,Description="parent">',
            '##INFO=<ID=submitter_variant_call_id,Number=.,Type=Integer,Description="submitter variant call id">',
            '##INFO=<ID=remap_score,Number=.,Type=Float,Description="remap score">',
            '##INFO=<ID=Variant_seq,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">',
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
            '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">',
            '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">',
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
            '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">',
            '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">'
        ]


    def test_get_alt(self):
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	81	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=77,78;End_range=80,81;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6], f_list[7], f_list[8])
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        info_attribute_dict = read_info_attributes(self.info_attribute_input_file)

        symbolic_allele_dictionary = read_sequence_ontology_symbolic_allele(self.symbolic_allele_file)
        assembly_file = self.assembly

        # standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        # merging
        standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }
        # Dictionary for all possible VCF meta-information lines
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT")
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO")
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER")
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT")

        all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        v = VcfLine(line_object,
                    info_attribute_dict,
                    symbolic_allele_dictionary,
                    assembly_file,
                    standard_lines_dictionary,
                    all_possible_lines_dictionary)
        alt_allele = v.get_alt(standard_lines_dictionary, all_possible_lines_dictionary)
        assert alt_allele == '<DEL>'

    def test_generate_vcf_metainformation(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)

        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          self.assembly)
        (unique_pragmas_to_add, sample_names,
         unique_alt_lines_to_add, unique_info_lines_to_add,
         unique_filter_lines_to_add, unique_format_lines_to_add) = generate_vcf_metainformation(
            gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary
        )
        print(unique_pragmas_to_add)
        assert unique_pragmas_to_add == ['##fileformat=VCFv4.4', '##gff-version=3', '##source=DGVa', '##gvf-version=1.06',
                                         '##species=http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##fileDate=2015-07-15',
                                         '##genome-build=NCBIGRCz10', '##Study_accession=nstd62', '##Study_type=Control Set',
                                         '##Display_name=Brown_et_al_2012', '##Publication_year=2012',
                                         '##Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants',
                                         '##Assembly_name=GRCz10', '##subject=subject_name=Wilds2-3', '##subject=subject_name=Zon9', '##subject=subject_name=JenMale7;subject_sex=Male',
                                         '##subject=subject_name=JenMale6;subject_sex=Male', '##sample=sample_name=JenMale6;subject_name=JenMale6', '##sample=sample_name=Wilds2-3;subject_name=Wilds2-3',
                                         '##sample=sample_name=Zon9;subject_name=Zon9', '##sample=sample_name=JenMale7;subject_name=JenMale7']


        assert unique_alt_lines_to_add == [
            '"##ALT=<ID=DEL,Description=""Deletion"">"',
            '"##ALT=<ID=DUP,Description=""Duplication"">"'
        ]
        print("unique_info_lines_to_add", unique_info_lines_to_add)
        assert unique_info_lines_to_add ==  [
            '##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier.">',
            '##INFO=<ID=Name,Number=.,Type=String,Description="name">',
            '##INFO=<ID=Alias,Number=.,Type=String,Description="A secondary name.">',
            '##INFO=<ID=variant_call_so_id,Number=.,Type=String,Description="variant call SO id">',
            '##INFO=<ID=parent,Number=.,Type=String,Description="parent">',
            '##INFO=<ID=submitter_variant_call_id,Number=.,Type=Integer,Description="submitter variant call id">',
            '##INFO=<ID=remap_score,Number=.,Type=Float,Description="remap score">',
            '##INFO=<ID=Variant_seq,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">',
            '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
            '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">',
            '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">',
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'
        ]

    def test_generate_vcf_header_line(self):
        header_fields = generate_vcf_header_line(['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'])
        assert header_fields == '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7'

    def test_populate_sample_formats(self):
        sample_name_format_value = populate_sample_formats(['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'])
        assert sample_name_format_value == {
            'JenMale6': 'sampleFORMAThere',
            'Wilds2-3': 'sampleFORMAThere',
            'Zon9': 'sampleFORMAThere',
            'JenMale7': 'sampleFORMAThere'
        }

    def test_format_sample_values(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        # standard structured meta-information lines for this VCF file
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list,
                                                                          self.assembly)
        unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add, unique_filter_lines_to_add, unique_format_lines_to_add = generate_vcf_metainformation(
            gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary)
        sample_name_format_value = populate_sample_formats(samples)
        sample_format_values_string = format_sample_values(sample_name_format_value)
        assert isinstance(sample_format_values_string, str)

    def test_format_vcf_datalines(self):
        gvf_pragmas, gvf_non_essential, gvf_lines_obj_list = read_in_gvf_file(self.input_file)
        header_standard_lines_dictionary, vcf_data_lines, list_of_vcf_objects = gvf_features_to_vcf_objects(gvf_lines_obj_list, self.assembly)
        (
            unique_pragmas_to_add, samples, unique_alt_lines_to_add, unique_info_lines_to_add,
            unique_filter_lines_to_add, unique_format_lines_to_add
         ) = generate_vcf_metainformation(gvf_pragmas, gvf_non_essential, list_of_vcf_objects, header_standard_lines_dictionary)
        formatted_vcf_datalines = format_vcf_datalines(list_of_vcf_objects, samples)
        assert formatted_vcf_datalines == [
            'chromosome1\t1\t1\tAC\t<DEL>\t.\t.\tEND=1;SVLEN=1\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t',
            'chromosome1\t76\t1\tTAA\t<DEL>\t.\t.\tEND=78;SVLEN=1;IMPRECISE;CIPOS=776537,776837;CIEND=776537,776837\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t',
            'chromosome1\t126\t12\tCGTACGGTACG\t<DEL>\t.\t.\tEND=131;SVLEN=5\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t',
            'chromosome1\t127\t13\tGTACGTACG\t<DUP>\t.\t.\tEND=131;SVLEN=4\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t',
            'chromosome1\t127\t14\tGTACGTACG\t<DUP>\t.\t.\tEND=131;SVLEN=4\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t',
            'chromosome1\t127\t14\tGTT\t<DUP>\t.\t.\tEND=128;SVLEN=1\tpending\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\tsampleFORMAThere\t'
        ]

    def test_generate_custom_unstructured_metainfomation_line(self):
        formatted_string = generate_custom_unstructured_meta_line("test_string_key", "test_string_value")
        assert formatted_string == "##test_string_key=test_string_value"

if __name__ == '__main__':
    unittest.main()
