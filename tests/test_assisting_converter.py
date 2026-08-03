import os
import unittest

from convert_gvf_to_vcf.gvf_attributes_converter import generate_custom_structured_meta_line, \
    GvfAttributeParser, GvfAttributeTransformer
from convert_gvf_to_vcf.lookup import Lookup
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestAssistingConverter(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.input_file = os.path.join(self.paths.test_dir, "input", "zebrafish.gvf")
        self.input_folder_parent = self.paths.package_dir
        # the inputs below are INFO attribute files
        self.etc_folder =  self.paths.etc_dir
        self.assembly = os.path.join(self.paths.test_dir,"input", "zebrafish.fa")
        self.output_file = os.path.join(self.paths.test_dir,"input", "a.vcf")
        self.reference_lookup = Lookup(self.assembly, self.paths)
        #
        self.mapping_dictionary = self.reference_lookup.mapping_attribute_dict
        self.field_lines_dictionary = {'ALT': [], 'INFO': [], 'FILTER': [], 'FORMAT': []}
        self.all_possible_lines_dictionary = {
            'ALT': {
                'DUP': '##ALT=<ID=DUP,Description="Duplication">',
                'DEL': '##ALT=<ID=DEL,Description="Deletion">'
            },
            'INFO': {
                'ID': '##INFO=<ID=ID,Number=.,Type=String,Description="A unique identifier">',
                'NAME': '##INFO=<ID=NAME,Number=.,Type=String,Description="Name">',
                'ALIAS': '##INFO=<ID=ALIAS,Number=.,Type=String,Description="Secondary Name">',
                'VARCALLSOID': '##INFO=<ID=VARCALLSOID,Number=.,Type=String,Description="Variant call Sequence ontology ID">',
                'SVCID': '##INFO=<ID=SVCID,Number=.,Type=Integer,Description="submitter variant call ID">',
                'VARSEQ': '##INFO=<ID=VARSEQ,Number=.,Type=String,Description="Alleles found in an individual (or group of individuals).">',
                'END': '##INFO=<ID=END,Number=1,Type=Integer,Description="End position on CHROM (used with symbolic alleles; see below) or End position of the longest variant described in this record">',
                'SVLEN': '##INFO=<ID=SVLEN,Number=A,Type=String,Description="Length of structural variant">',
                'IMPRECISE': '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
                'CIPOS': '##INFO=<ID=CIPOS,Number=.,Type=Integer,Description="Confidence interval around POS for symbolic structural variants">',
                'CIEND': '##INFO=<ID=CIEND,Number=.,Type=Integer,Description="Confidence interval around END for symbolic structural variants">',
                'AC': '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">',
                'DBXREF': '##INFO=<ID=DBXREF,Number=.,Type=String,Description="A database cross-reference">',
                'AD': '##INFO=<ID=AD,Number=R,Type=Integer,Description="Total read depth for each allele">'
            },
            'FILTER': {
            },
            'FORMAT': {
                'AD': '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
                'GT': '##FORMAT=<ID=GT,Number=1,Type=String,Description="Zygosity">'
            }
        }

    def test_generate_custom_structured_meta_line(self):
        field = "INFO"
        idkey = "END"
        number = "1"
        data_type = "Integer"
        description = "Testing the end position"
        # Testing None
        optional_data_none = None
        output_for_none = generate_custom_structured_meta_line(
            field,
            idkey,
            number,
            data_type,
            description,
            optional_data_none
        )
        expected_output_for_none = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Testing the end position\">"
        assert output_for_none == expected_output_for_none
        # Testing optional data
        optional_data = {"test1": "A"}
        output_for_optional = generate_custom_structured_meta_line(
            field,
            idkey,
            number,
            data_type,
            description,
            optional_data
        )
        expected_output_for_optional = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Testing the end position\",test1=\"A\">"
        assert output_for_optional == expected_output_for_optional

    def test_get_gvf_attributes(self):
        col9 = "ID=1;Name=nssv1412199;Alias=CNV28955;parent=nsv811094;Start_range=.,1;End_range=2,.;sample_name=Wilds2-3"
        gvf_parser = GvfAttributeParser(col9)
        expected_dictionary = {'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'parent': 'nsv811094', 'Start_range': ['.', '1'], 'End_range': ['2', '.'], 'sample_name': 'Wilds2-3'}
        assert gvf_parser._attributes == expected_dictionary

    def test_convert_gvf_attributes_to_vcf_values(self):
        col9 = "ID=1;Name=nssv1412199;Alias=CNV28955;parent=nsv811094;Start_range=.,1;End_range=2,.;sample_name=Wilds2-3;Genotype=0:1"

        gvf_parser = GvfAttributeParser(col9)

        transformer = GvfAttributeTransformer(
            mapping_attribute_dict=self.mapping_dictionary,
            field_lines_dictionary=self.field_lines_dictionary,
            all_possible_lines_dictionary=self.all_possible_lines_dictionary
        )
        gvf_attribute_dictionary, vcf_info_values, vcf_format_values = transformer.convert_gvf_attributes_to_vcf_values(
            gvf_parser)


        # testing gvf_attribute_dictionary
        expected_gvf_attribute_dictionary = {'ID': '1', 'Name': 'nssv1412199', 'Alias': 'CNV28955', 'parent': 'nsv811094', 'Start_range': ['.', '1'], 'End_range': ['2', '.'], 'sample_name': 'Wilds2-3', 'Genotype': '0:1'}
        assert gvf_attribute_dictionary == expected_gvf_attribute_dictionary
        # testing vcf_info_values
        expected_info_values = {'ID': '1', 'NAME': 'nssv1412199', 'ALIAS': 'CNV28955', 'PARENT': 'nsv811094', 'SAMPLENAME': 'Wilds2-3'}
        assert vcf_info_values == expected_info_values
        # testing vcf_format_values
        expected_format_values = {'Wilds2-3': {'GT': '0:1'}}
        assert vcf_format_values == expected_format_values

    def test_apply_genotype_inference(self):
        vcf_info_values = {}
        vcf_format_values = {}

        gvf_parser = GvfAttributeParser("sample_name=Wilds2-3")
        transformer = GvfAttributeTransformer(None, None, None)
        transformer.apply_genotype_inference(gvf_parser, vcf_format_values, vcf_info_values)

        expected_format_values = {'Wilds2-3': {'GT': '1/.'}}
        assert vcf_format_values == expected_format_values

    def test_process_vcf_fields(self):
        gvf_parser = GvfAttributeParser("ID=1;sample_name=Wilds1;Genotype=0:1")
        vcf_info_values = {}
        vcf_format_values = {}

        info_field_values = {"INFO": {"FieldKey": "ID", "Number": "1", "Type": "String", "Description": "ID"}}
        transformer = GvfAttributeTransformer(self.mapping_dictionary, self.field_lines_dictionary, self.all_possible_lines_dictionary)
        transformer.process_vcf_fields("ID", "INFO", info_field_values, gvf_parser, vcf_format_values, vcf_info_values)
        assert vcf_info_values == {"ID": "1"}

        format_field_values = {"FORMAT": {"FieldKey": "GT"}}
        transformer.process_vcf_fields("Genotype", "FORMAT", format_field_values, gvf_parser, vcf_format_values, vcf_info_values)
        assert vcf_format_values == {"Wilds1": {"GT": "0:1"}}

    def test_infer_genotype(self):
        transformer = GvfAttributeTransformer(self.mapping_dictionary, self.field_lines_dictionary, self.all_possible_lines_dictionary)

        gvf_attributes = {"sample_name": "Wilds2-3"}
        vcf_format_values = {}

        transformer.infer_genotype(gvf_attributes, vcf_format_values)

        expected_format_values = {'Wilds2-3': {'GT': '1/.'}}
        assert vcf_format_values == expected_format_values

    def test_determine_vcf_type(self):
        transformer = GvfAttributeTransformer(self.mapping_dictionary, self.field_lines_dictionary, self.all_possible_lines_dictionary)

        vcf_info_values = {}
        vcf_format_values = {'Wilds2-3': {'GT': '0/1'}}
        result = transformer.determine_vcf_type(vcf_info_values, vcf_format_values)
        assert result == "GENOTYPED"

        vcf_info_values = {'AF': '0.5', 'AC': '1', 'AN': '2'}
        vcf_format_values = {}
        result = transformer.determine_vcf_type(vcf_info_values, vcf_format_values)
        assert result == "POPULATED"

        vcf_info_values = {}
        vcf_format_values = {}
        result = transformer.determine_vcf_type(vcf_info_values, vcf_format_values)
        assert result == "SITES-ONLY"

    def test_process_vcf_format_field(self):

        transformer = GvfAttributeTransformer(self.mapping_dictionary, self.field_lines_dictionary, self.all_possible_lines_dictionary)

        attrib_key = "Genotype"
        field = "FORMAT"
        field_values = {"FORMAT": {"FieldKey": "GT"}}
        gvf_attribute_dictionary = {"sample_name": "Wilds2-3", "Genotype": "0:1"}
        vcf_format_values = {}

        transformer.process_vcf_format_field(attrib_key, field, field_values,
                                             gvf_attribute_dictionary, vcf_format_values)


        expected_header = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Zygosity">'
        assert self.field_lines_dictionary['FORMAT'] == [expected_header]


        expected_format_values = {'Wilds2-3': {'GT': '0:1'}}
        assert vcf_format_values == expected_format_values


    def test_process_vcf_info_field(self):
        transformer = GvfAttributeTransformer(self.mapping_dictionary, self.field_lines_dictionary, self.all_possible_lines_dictionary)

        attrib_key = "ID"
        field = "INFO"
        field_values = {
            "INFO": {
                "FieldKey": "ID",
                "Number": ".",
                "Type": "String",
                "Description": "A unique identifier"
            }
        }
        gvf_attribute_dictionary = {"ID": "1"}
        vcf_info_values = {}

        transformer.process_vcf_info_field(attrib_key, field, field_values,
                                           gvf_attribute_dictionary, vcf_info_values)


        assert len(self.field_lines_dictionary['INFO']) == 1
        assert "ID=ID" in self.field_lines_dictionary['INFO'][0]

        assert vcf_info_values == {"ID": "1"}

