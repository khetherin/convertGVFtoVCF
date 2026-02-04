import unittest
import os

from convert_gvf_to_vcf.conversionstatistics import FileStatistics, StatisticsPayload
from convert_gvf_to_vcf.convertGVFtoVCF import convert_gvf_pragmas_for_vcf_header, generate_vcf_header_line
from convert_gvf_to_vcf.lookup import Lookup


class TestConversionStatistics(unittest.TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.gvf_file = os.path.join(input_folder, "input", "zebrafish.gvf")
        self.vcf_file = os.path.join(input_folder, "input", "a.vcf")
        self.assembly = os.path.join(input_folder, "input", "zebrafish.fa")
        self.reference_lookup = Lookup(self.assembly)
        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06', '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955', '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        gvf_pragma_comments = ['#Study_accession: nstd62', '#Study_type: Control Set', '#Display_name: Brown_et_al_2012', '#Publication: PMID=22203992;Journal=Proceedings of the National Academy of Sciences of the United States of America;Paper_title=Extensive genetic diversity and substructuring among zebrafish strains revealed through copy number variant analysis.;Publication_year=2012', '#Study: First_author=Kim Brown;Description=Comparative genomic hybridization analysis of 3 laboratory and one wild zebrafish populations for Copy Number Variants', '#Assembly_name: GRCz10', '#subject: subject_name=Wilds2-3', '#subject: subject_name=Zon9', '#subject: subject_name=JenMale7;subject_sex=Male', '#subject: subject_name=JenMale6;subject_sex=Male', '#sample: sample_name=JenMale6;subject_name=JenMale6', '#sample: sample_name=Wilds2-3;subject_name=Wilds2-3', '#sample: sample_name=Zon9;subject_name=Zon9', '#sample: sample_name=JenMale7;subject_name=JenMale7', '#testing_unknown_pragma']
        _, unique_sample_name =convert_gvf_pragmas_for_vcf_header(gvf_pragma, gvf_pragma_comments, self.reference_lookup)
        self.samples = unique_sample_name
        self.payload = StatisticsPayload(samples=['JenMale6', 'Wilds2-3', 'Zon9', 'JenMale7'],
                                         vcf_header_fields='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tJenMale6\tWilds2-3\tZon9\tJenMale7',
                                         gvf_line_count=7, vcf_line_count=5, vcf_merge_count=2,
                                         gvf_chromosome_count={'chromosome1': 7},
                                         vcf_chromosome_count={'chromosome1': 5},
                                         gvf_sv_so_term_count={'copy_number_loss': 3, 'copy_number_gain': 4},
                                         vcf_alt_count={'<DEL>': 3, '<DUP>': 2}, vcf_alt_missing=0,
                                         vcf_info_count = {'ID': 3, 'NAME': 5, 'ALIAS': 5, 'VARCALLSOID': 5, 'PARENT': 5, 'SVCID': 5, 'SAMPLENAME': 5,
                            'REMAP': 5, 'VARSEQ': 5, 'END': 5, 'IMPRECISE': 5, 'SVLEN': 5, 'SVCLAIM': 5, 'CIPOS': 4,
                            'CIEND': 4, 'AC': 2, 'AD': 1, 'DBXREF': 1}
                                         )

    def test_get_info_count(self):
        statistics_summariser = FileStatistics(self.vcf_file, self.payload)
        # print(len(statistics_summariser.get_info_count(self.payload)))
        assert len(statistics_summariser.get_info_count(self.payload)) == 18
        assert statistics_summariser.get_info_count(self.payload)["ID"] == 3
        assert statistics_summariser.get_info_count(self.payload)["NAME"] == 5



    def test_get_alt_missing(self):
        statistics_summariser = FileStatistics(self.vcf_file, self.payload)
        assert statistics_summariser.get_alt_alleles_missing(self.payload) == 0

    def test_get_variant_lines(self):
        #GVF
        statistics_summariser = FileStatistics(self.gvf_file, self.payload)
        assert statistics_summariser.get_variant_lines(self.payload) == 7
        # VCF
        statistics_summariser = FileStatistics(self.vcf_file, self.payload)
        assert statistics_summariser.get_variant_lines(self.payload) == 5


    def test_get_sample_number_missing(self):
        # VCF
        is_missing_format_value = True # all format fields are missing
        header_fields = generate_vcf_header_line(is_missing_format=is_missing_format_value, samples=self.samples)
        payload = StatisticsPayload(
            samples=self.samples,
            vcf_header_fields=header_fields,
            gvf_line_count=None,
            vcf_line_count=None,
            vcf_merge_count=None,
            gvf_chromosome_count=None,
            vcf_chromosome_count=None,
            gvf_sv_so_term_count= None,
            vcf_alt_count=None,
            vcf_alt_missing=None,
            vcf_info_count=None
        )
        statistics_summariser = FileStatistics(self.vcf_file, payload)
        assert statistics_summariser.get_sample_number_missing(self.samples, header_fields) == 4


    def test_get_sample_number(self):

        # GVF: set up to get the samples
        statistics_summariser = FileStatistics(self.gvf_file, self.payload)
        sample_number_gvf = statistics_summariser.get_sample_number(self.samples, vcf_header_fields=self.payload.vcf_header_fields)
        assert sample_number_gvf == 4
        # VCF
        is_missing_format_value = False # all format fields are present
        # is_missing_format_value = True # all format fields are missing
        header_fields = generate_vcf_header_line(is_missing_format=is_missing_format_value, samples=self.samples)
        payload = StatisticsPayload(
            samples=self.samples,
            vcf_header_fields=header_fields,
            gvf_line_count=None,
            vcf_line_count=None,
            vcf_merge_count=None,
            gvf_chromosome_count = None,
            vcf_chromosome_count=None,
            gvf_sv_so_term_count=None,
            vcf_alt_count=None,
            vcf_alt_missing=None,
            vcf_info_count=None
        )
        statistics_summariser = FileStatistics(self.vcf_file, payload)
        assert statistics_summariser.get_sample_number(self.samples, header_fields) == 4

    def test_get_file_md5(self):
        statistics_summariser = FileStatistics(self.gvf_file, self.payload)
        md5 = statistics_summariser.get_file_md5()
        assert md5 == "cd842fc345bd990d4b14fd80c1d08f8e"

    def test_get_file_version(self):
        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06',
                      '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955',
                      '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        gvf_pragma_comments = ['#Study_accession: nstd62', '#Study_type: Control Set']

        unique_converted_pragmas, unique_sample_name = convert_gvf_pragmas_for_vcf_header(gvf_pragma,
                                                                                          gvf_pragma_comments,
                                                                                          self.reference_lookup)
        statistics_summariser = FileStatistics(self.gvf_file, self.payload)
        version_gvf = statistics_summariser.get_file_version(gvf_pragma)
        assert version_gvf == "##gvf-version 1.06"

        # testing for VCF
        is_missing_format_value = False # all format fields are present
        header_fields = generate_vcf_header_line(is_missing_format=is_missing_format_value, samples=self.samples)
        payload = StatisticsPayload(
            samples=self.samples,
            vcf_header_fields=header_fields,
            gvf_line_count=None,
            vcf_line_count=None,
            vcf_merge_count=None,
            gvf_chromosome_count=None,
            vcf_chromosome_count=None,
            gvf_sv_so_term_count=None,
            vcf_alt_count=None,
            vcf_alt_missing=None,
            vcf_info_count=None
        )
        statistics_summariser = FileStatistics(self.vcf_file, payload)
        version_vcf = statistics_summariser.get_file_version(unique_converted_pragmas)
        assert version_vcf == "##fileformat=VCFv4.4"


    def test___str__(self):
        statistics_summariser = FileStatistics(self.gvf_file, self.payload)

        gvf_pragma = ['##gff-version 3', '##gvf-version 1.06',
                      '##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7955',
                      '##file-date 2015-07-15', '##genome-build NCBI GRCz10']
        statistics_summariser.get_file_version(gvf_pragma)
        statistics_summary = statistics_summariser.__str__()

        (file_name_string,
         file_path,
         file_extension,
         file_size_string,
         timestamp_string,
         file_version_string,
         file_md5_string,
         end_string,
         sample_number,
         gvf_chromsome_count,
         vcf_chromosome_count,
         gvf_sv_so_count,
         _
         ) = statistics_summary.split("\n")
        assert file_name_string == "File name = zebrafish.gvf"
        assert file_size_string == "File size = 2971 bytes"
        assert timestamp_string == "Last mod Timestamp = 2026-01-27 11:32:30.488227"
        assert file_version_string == "Version = ##gvf-version 1.06"
        assert file_md5_string == "md5 string = cd842fc345bd990d4b14fd80c1d08f8e"
        assert end_string == "====="
        print(gvf_sv_so_count)

if __name__ == '__main__':
    unittest.main()

