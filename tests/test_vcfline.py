import os
import unittest

from Bio.PDB import refmat

from convert_gvf_to_vcf.convertGVFtoVCF import generate_vcf_header_structured_lines
from convert_gvf_to_vcf.gvffeature import GvfFeatureline
from convert_gvf_to_vcf.vcfline import VcfLine, VcfLineBuilder
from convert_gvf_to_vcf.lookup import Lookup

class TestVcfLineBuilder(unittest.TestCase):
    def setUp(self):
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'convert_gvf_to_vcf'))
        self.etc_folder = os.path.join(self.input_folder_parent, "etc")
        input_folder = os.path.dirname(__file__)
        # Prepare Inputs
        self.input_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.input_folder_parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'convert_gvf_to_vcf'))
        # Prepare Outputs
        self.output_file = os.path.join(input_folder, "input", "a.vcf")
        # Prepare References
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly)
        # Set up GVF line object
        gvf_feature_line = "chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."
        f_list = gvf_feature_line.split("\t")
        gvf_line_object = GvfFeatureline(f_list[0], f_list[1], f_list[2], f_list[3], f_list[4], f_list[5], f_list[6],
                                         f_list[7], f_list[8])
        # Set up of data structures
        # Dictionary of standard structured meta-information lines for this VCF file
        lines_standard_alt = []
        lines_standard_info = []
        lines_standard_filter = []
        lines_standard_format = []
        self.standard_lines_dictionary = {
            "ALT": lines_standard_alt,
            "INFO": lines_standard_info,
            "FILTER": lines_standard_filter,
            "FORMAT": lines_standard_format,
        }

        # Dictionary for all possible VCF meta-information lines
        all_possible_info_lines = generate_vcf_header_structured_lines("INFO",
                                                                       self.reference_lookup.mapping_attribute_dict)
        all_possible_alt_lines = generate_vcf_header_structured_lines("ALT",
                                                                      self.reference_lookup.mapping_attribute_dict)
        all_possible_filter_lines = generate_vcf_header_structured_lines("FILTER",
                                                                         self.reference_lookup.mapping_attribute_dict)
        all_possible_format_lines = generate_vcf_header_structured_lines("FORMAT",
                                                                         self.reference_lookup.mapping_attribute_dict)
        self.all_possible_lines_dictionary = {
            "ALT": all_possible_alt_lines,
            "INFO": all_possible_info_lines,
            "FILTER": all_possible_filter_lines,
            "FORMAT": all_possible_format_lines,
        }
        ordered_list_of_samples = ['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7']
        self.vcf_builder = VcfLineBuilder(self.standard_lines_dictionary, self.all_possible_lines_dictionary,
                                          self.reference_lookup, ordered_list_of_samples)

    def test_build_vcf_line(self):
        # Set up GVF line object
        gvf_attributes = 'ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=.'
        gvf_line_object = GvfFeatureline('chromosome1', 'DGVa', 'copy_number_loss', '77', '78', '.', '+', '.', gvf_attributes )
        vcf_line = self.vcf_builder.build_vcf_line(gvf_line_object)
        assert vcf_line.chrom == 'chromosome1'
        assert vcf_line.pos == 76
        assert vcf_line.info_dict['NAME'] == 'nssv1412199'
        assert vcf_line.info_dict['ALIAS'] == 'CNV28955'

    def test_add_padded_base(self):
        padded_base, pos, ref, alt = self.vcf_builder.add_padded_base(chrom='chromosome1', pos=79, end=80, ref="A", alt='.', placed_before=True)
        assert padded_base == 'C'
        assert pos  == 78
        assert ref == 'CA'
        assert alt  == '.'
        padded_base, pos, ref, alt = self.vcf_builder.add_padded_base(chrom='chromosome1', pos=1, end=2, ref="A", alt='.', placed_before=False)
        assert padded_base == 'C'
        assert pos  == 1
        assert ref == 'AC'
        assert alt  == '.'

    def test_get_ref(self):
        reference_allele = self.vcf_builder.get_ref(vcf_value_from_gvf_attribute={}, chrom='chromosome1', pos=1, end=2)
        assert len(reference_allele) != 0
        assert reference_allele == 'A'

    def test_create_coordinate_range(self):
        # with range
        vcf_value_from_gvf_attribute_with_range={'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'variant_call_so_id': 'SO:0001743',
         'parent': 'nsv811094', 'Start_range': ['.', '1'], 'End_range': ['2', '.'],
         'submitter_variant_call_id': 'CNV28955', 'sample_name': 'Wilds2-3', 'remap_score': '.98857',
         'Variant_seq': '.'}
        start_range_lower_bound, start_range_upper_bound, end_range_lower_bound, end_range_upper_bound = self.vcf_builder.create_coordinate_range(vcf_value_from_gvf_attribute_with_range, pos=int(1),end=int(2))
        assert start_range_lower_bound == "."
        assert start_range_upper_bound == "1"
        assert end_range_lower_bound == "2"
        assert end_range_upper_bound == "."
        # without range
        vcf_value_from_gvf_attribute_without_range={'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'variant_call_so_id': 'SO:0001743',
         'parent': 'nsv811094',
         'submitter_variant_call_id': 'CNV28955', 'sample_name': 'Wilds2-3', 'remap_score': '.98857',
         'Variant_seq': '.'}
        start_range_lower_bound, start_range_upper_bound, end_range_lower_bound, end_range_upper_bound = self.vcf_builder.create_coordinate_range(
            vcf_value_from_gvf_attribute_with_range, pos=int(1), end=int(2))
        assert start_range_lower_bound == "."
        assert start_range_upper_bound == "1"
        assert end_range_lower_bound == "2"
        assert end_range_upper_bound == "."

    def test_generate_symbolic_allele(self):
        # imprecise variant
        vcf_value_from_gvf_attribute = {"Variant_seq":".", "Start_range":".,776614","End_range":"786127,."}
        symbolic_allele, info_dict = self.vcf_builder.generate_symbolic_allele(vcf_value_from_gvf_attribute, pos=76, end=77, length=1, ref='.', so_type='copy_number_loss')
        assert symbolic_allele == '<DEL>'
        assert info_dict == {'END': '77', 'IMPRECISE': 'IMPRECISE', 'CIPOS': '.,776538', 'CIEND': '786050,.', 'SVLEN': '1'}
        # precise variant
        vcf_value_from_gvf_attribute = {"Variant_seq": "."}
        symbolic_allele, info_dict = self.vcf_builder.generate_symbolic_allele(vcf_value_from_gvf_attribute, pos=76,
                                                                               end=77, length=1, ref='.',
                                                                               so_type='copy_number_loss')
        assert symbolic_allele == '<DEL>'
        assert info_dict == {'END': '76', 'IMPRECISE': None, 'CIPOS': None, 'CIEND': None, 'SVLEN': '1'}


    def test_generate_info_field_symbolic_allele(self):
        end = 100
        end_range_lower_bound = 100
        end_range_upper_bound = 200
        info_end_value = None
        length = 13
        pos = 1
        ref = "A" * 13
        start_range_lower_bound = 1
        start_range_upper_bound = 2
        symbolic_allele = "<INS>"
        info_dict, is_imprecise = self.vcf_builder.generate_info_field_symbolic_allele(
            end, end_range_lower_bound, end_range_upper_bound,
            length, pos, ref,
            start_range_lower_bound, start_range_upper_bound,
            symbolic_allele
        )
        assert info_dict == {'END': '13', 'IMPRECISE': 'IMPRECISE', 'CIPOS': '0,1', 'CIEND': '0,100', 'SVLEN': '13'}
        assert is_imprecise is True

    def test_generate_info_field_for_imprecise_variant(self):
        end = 100
        end_range_lower_bound = 100
        end_range_upper_bound = 200
        info_end_value = None
        length = 13
        pos = 1
        ref = "A"*13
        start_range_lower_bound = 1
        start_range_upper_bound = 2
        symbolic_allele = "<INS>"
        info_ciend_value, info_cipos_value, info_end_value, info_imprecise_value, is_imprecise = self.vcf_builder.generate_info_field_for_imprecise_variant(
            end, end_range_lower_bound, end_range_upper_bound,
            info_end_value,
            length, pos, ref,
            start_range_lower_bound, start_range_upper_bound,
            symbolic_allele
        )
        assert info_ciend_value == "0,100"
        assert info_cipos_value == "0,1"
        assert info_end_value == "13"
        assert info_imprecise_value == "IMPRECISE"
        assert is_imprecise is True

    def test_generate_info_field_for_precise_variant(self):
        pos = 34
        ref = "AAATTTGGG"
        info_end_value, is_imprecise = self.vcf_builder.generate_info_field_for_precise_variant(pos, ref)
        assert info_end_value == "42"
        assert is_imprecise is False

    def test_calculate_CIPOS_and_CIEND(self):
        # known imprecise variant
        end=133
        end_range_lower_bound=122
        end_range_upper_bound=133
        pos=1
        start_range_lower_bound=1
        start_range_upper_bound=10
        ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound = self.vcf_builder.calculate_CIPOS_and_CIEND(end, end_range_lower_bound, end_range_upper_bound, pos, start_range_lower_bound, start_range_upper_bound)
        assert (ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound) == (-11, 0, 0, 9)
        # unknown precise variant
        end=122
        end_range_lower_bound=122
        end_range_upper_bound="."
        pos=1
        start_range_lower_bound="."
        start_range_upper_bound=10
        ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound = self.vcf_builder.calculate_CIPOS_and_CIEND(end, end_range_lower_bound, end_range_upper_bound, pos, start_range_lower_bound, start_range_upper_bound)
        assert (ciend_lower_bound, ciend_upper_bound, cipos_lower_bound, cipos_upper_bound) == (0, ".", ".", 9)

    def test_get_alt(self):
            'chromosome1	DGVa	copy_number_loss	77	78	.	+	.	ID=1;Name=nssv1412199;Alias=CNV28955;variant_call_so_id=SO:0001743;parent=nsv811094;Start_range=.,776614;End_range=786127,.;submitter_variant_call_id=CNV28955;sample_name=Wilds2-3;remap_score=.98857;Variant_seq=."'
            vcf_value_from_gvf_attribute = {"Variant_seq": ".", "Start_range": ".,776614", "End_range": "786127,."}
            pos, ref, alt, info_dict = self.vcf_builder.get_alt(vcf_value_from_gvf_attribute, chrom='chromosome1',
                                                                pos=77, end=78, length=1, ref='',
                                                                so_type='copy_number_loss')
            assert pos == 76
            assert ref == 'T'
            assert alt == '<DEL>'
            assert info_dict == {'END': '76', 'IMPRECISE': None, 'CIPOS': None, 'CIEND': None, 'SVLEN': '1'}

    def test_check_ref(self):
        assert self.vcf_builder.check_ref('A') == 'A'
        assert self.vcf_builder.check_ref('T') == 'T'
        assert self.vcf_builder.check_ref('C') == 'C'
        assert self.vcf_builder.check_ref('G') == 'G'
        assert self.vcf_builder.check_ref('TGCR') == 'TGCA'

class TestVcfline(unittest.TestCase):
    def setUp(self):
        self.vcf_line = VcfLine(chrom='Chromosome1', pos='76', id='1', ref='T', alt='<DEL>', qual='.', filter='.',
                                info_dict={'NAME': 'nssv1412199', 'SVLEN': '1'},
                                vcf_values_for_format={'sample1':{'GT':'0/1'}}, order_sample_names=['sample1'])
        self.other_vcf_line = VcfLine(chrom='Chromosome1', pos='76', id='1', ref='T', alt='<DEL>', qual='.', filter='.',
                                info_dict={'NAME': 'nssv9912199', 'SVLEN': '4'},
                                vcf_values_for_format={'sample2':{'GT':'0/1'}}, order_sample_names=['sample2'])

    def test__str__(self):
        assert str(self.vcf_line) == 'Chromosome1\t76\t1\tT\t<DEL>\t.\t.\tNAME=nssv1412199;SVLEN=1\tGT\t0/1'

    def test__eq__(self):
        vcf_line1 = VcfLine(
            chrom='Chromosome1', pos='76', id='1', ref='T', alt='<DEL>', qual='.', filter='.',
            info_dict={},vcf_values_for_format={}
        )
        vcf_line2 = VcfLine(
            chrom='Chromosome1', pos='76', id='2', ref='T', alt='', qual='.', filter='.',
            info_dict={}, vcf_values_for_format={}
        )
        vcf_line3 = VcfLine(
            chrom='Chromosome1', pos='77', id='4', ref='A', alt='ATT', qual='.', filter='.',
            info_dict={}, vcf_values_for_format={}
        )
        assert vcf_line1 == vcf_line2
        assert vcf_line1 != vcf_line3

    def test_merge_and_add(self):
        # testing merge for different elements
        merged_string = self.vcf_line.merge_and_add('1', '2', ';')
        assert merged_string == '1;2'
        # testing non merge for same elements
        non_merged_string = self.vcf_line.merge_and_add("1", "1", ";")
        assert non_merged_string == "1"

    def test_order_format_keys(self):
        assert self.vcf_line._order_format_keys(['DP', 'PL', 'GT', 'GQ']) == ['GT', 'DP', 'GQ', 'PL']

    def test_combine_format_values_by_sample_as_str(self):
        assert self.vcf_line.combine_format_values_by_sample_as_str() == '0/1'
        format_values = {
            'sample1':{'GT':'0/1', 'DP':5, 'PL':'50,0,12', 'GQ':'25'},
            'sample2': {'GT': '1/1', 'DP': '15', 'PL': '120,45,0', 'GQ': '50'},
            'sample3': {'GT': '0/0', 'DP': 2, 'PL': '0,10,10', 'GQ': '2'},
        }
        vcf_line1 = VcfLine(chrom='Chromosome1', pos='76', id='1', ref='T', alt='<DEL>', qual='.', filter='.',
                            info_dict={},
                            vcf_values_for_format=format_values,
                            order_sample_names=['sample1'])
        # only print sample1
        assert vcf_line1.combine_format_values_by_sample_as_str() == '0/1:5:25:50,0,12'
        vcf_line1.order_sample_names = ['sample1', 'sample2', 'sample3']
        # Now prints all 3 samples
        assert vcf_line1.combine_format_values_by_sample_as_str()  == '0/1:5:25:50,0,12\t1/1:15:50:120,45,0\t0/0:2:2:0,10,10'

    def test_fill_merge_dicts(self):
        merged_info_dict = {'END': '128', 'AD': '3', 'VARSEQ': '.', 'ALIAS': 'CNV5711', 'ID': '14', 'VARCALLSOID': 'SO:0001742', 'SVCID': 'CNV5711', 'DBXREF': 'mydata', 'NAME': 'nssv1388955', 'AC': '3', 'SVLEN': '1'}
        key = "VARCALLSOID"
        previous_line_info_value = "SO:0001742"
        current_line_info_value = "SO:1234567"
        merged_info_dict= VcfLine.fill_merge_dicts(self, merged_info_dict, key, previous_line_info_value, current_line_info_value)
        expected_merged_info_dict = {'END': '128', 'AD': '3', 'VARSEQ': '.', 'ALIAS': 'CNV5711', 'ID': '14', 'VARCALLSOID': 'SO:0001742,SO:1234567', 'SVCID': 'CNV5711', 'DBXREF': 'mydata', 'NAME': 'nssv1388955', 'AC': '3', 'SVLEN': '1'}
        assert merged_info_dict == expected_merged_info_dict

    def test_merge_info_dicts(self):
        VcfLine.merge_info_dicts(
            self.vcf_line, self.other_vcf_line
        )
        expected_merge_info_dict = {'NAME': 'nssv1412199,nssv9912199', 'SVLEN': '1,4'}
        assert self.vcf_line.info_dict == expected_merge_info_dict
        assert self.other_vcf_line.info_dict == expected_merge_info_dict

    def test_merge_vcf_values_for_format(self):
        self.vcf_line.merge_vcf_values_for_format(self.other_vcf_line)
        expected_merge_vcf_values_for_format = {'sample1': {'GT': '0/1'}, 'sample2': {'GT': '0/1'}}
        assert self.vcf_line.vcf_values_for_format == expected_merge_vcf_values_for_format

    def test_format_info_string(self):
        formatted_info_string = VcfLine.format_info_string(self.vcf_line)
        expected_info_string = "NAME=nssv1412199;SVLEN=1"
        assert formatted_info_string == expected_info_string

    def test_merge(self):
        list_of_samples = ['sample1']
        VcfLine.merge(self.vcf_line, self.other_vcf_line, list_of_samples)
        assert self.vcf_line.id == self.other_vcf_line.id
        assert self.vcf_line.alt == self.other_vcf_line.alt
        assert self.vcf_line.filter == self.other_vcf_line.filter
        assert self.vcf_line.info_dict == self.other_vcf_line.info_dict
        expected_union = {'sample1': {'GT': '0/1'}, 'sample2': {'GT': '0/1'}}
        vcfline_union = self.vcf_line.vcf_values_for_format | self.other_vcf_line.vcf_values_for_format
        assert vcfline_union == expected_union

    def test_keep(self):
        list_of_samples = ['sample1']
        VcfLine.keep(self, list_of_samples)
        assert self.vcf_line.order_sample_names == list_of_samples
