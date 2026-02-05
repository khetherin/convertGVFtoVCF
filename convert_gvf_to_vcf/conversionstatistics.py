"""
The purpose of this file is to calculate the (conversion) statistics of the GVF and VCF file.
"""
import os
from dataclasses import dataclass
from datetime import datetime
import hashlib
from typing import Optional, Union
from collections import Counter

# @dataclass
# class StatisticsPayload:
#     """ Data transfer object for statistics input"""
#     samples: list[str]
#     vcf_header_fields: Optional[Union[str, None]]
#     gvf_line_count: Optional[Union[int, None]]
#     vcf_line_count: Optional[Union[int, None]]
#     vcf_merge_count: Optional[Union[int, None]]
#     gvf_chromosome_count: Optional[Union[dict, None]]
#     vcf_chromosome_count: Optional[Union[dict, None]]
#     gvf_sv_so_term_count: Optional[Union[dict, None]]
#     vcf_alt_count: Optional[Union[dict, None]]
#     vcf_alt_missing: Optional[Union[int, None]]
#     vcf_info_count: Optional[Union[dict, None]]

class FileStatistics:
    """
    The responsibility of this class is to determine the statistics of a GVF or VCF file.
    """
    def __init__(self, gvf_file_path, gvf_pragmas, samples):
        self.gvf_file_path = gvf_file_path
        # file name
        self.gvf_file_name = os.path.basename(self.gvf_file_path)
        self.gvf_extension = os.path.splitext(self.gvf_file_path)[1]
        # file metadata
        self.gvf_metadata_stats = os.stat(self.gvf_file_path)
        self.gvf_file_size = self.gvf_metadata_stats.st_size
        gvf_last_modified_timestamp_raw = self.gvf_metadata_stats.st_mtime
        self.gvf_last_modified_timestamp = datetime.fromtimestamp(gvf_last_modified_timestamp_raw)
        self.gvf_version = self.find_version("##gvf-version", gvf_pragmas)
        self.gvf_md5 = self.get_file_md5(self.gvf_file_path)
        # biological samples
        self.sample_number = len(samples)
        # self.sample_number = self.get_sample_number(payload.samples, payload.vcf_header_fields)
        # self.sample_number_missing = self.get_sample_number_missing(payload.samples, payload.vcf_header_fields)
        # # chromsomes
        self.gvf_chromosome_count = Counter()
        self.vcf_chromosome_count = Counter()
        # # variants
        self.gvf_feature_line_count = 0
        self.vcf_data_line_count = 0
        # self.variant_lines = self.get_variant_lines(payload)
        self.vcf_number_of_merges = 0
        self.vcf_alt_alleles_count = Counter()
        # # TODO: number of additional ALT alleles
        # self.vcf_alt_missing = self.get_alt_alleles_missing(payload)
        # # TODO: number of SV types
        # self.sv_types = self.get_sv_types(payload)
        self.gvf_sv_so_term_count = Counter()
        # # TODO: number of precise variants
        # # TODO: number of imprecise variants
        # # attribute mapping
        # # TODO: observed GVF attributes (name and count)
        # # TODO: FORMAT tags (name and count)
        self.vcf_info_counter = Counter()
        # self.info_count = self.get_info_count(payload)
        # # TODO: dropped GVF attributes (assistingconverter.py:convert_gvf_attributes_to_vcf_values)
        # #convert > build_vcf_line> convert_gvf_attributes_to_vcf_values
        # # TODO: number of SO terms

    # def get_info_count(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.info_count = None
    #     if self.extension.endswith(".vcf"):
    #         self.info_count = payload.vcf_info_count
    #     return self.info_count
    #
    # def get_alt_alleles_missing(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.vcf_alt_missing = None
    #     if self.extension.endswith(".vcf"):
    #         self.vcf_alt_missing = payload.vcf_alt_missing
    #     return self.vcf_alt_missing
    #
    # def get_alt_alleles(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.vcf_alt_alleles = None;
    #     if self.extension.endswith(".vcf"):
    #         self.vcf_alt_alleles = payload.vcf_alt_count
    #     return self.vcf_alt_alleles
    #
    # def get_sv_types(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.sv_types = payload.gvf_sv_so_term_count
    #     if self.extension.endswith(".vcf"):
    #         #TODO: add the SO terms
    #         self.sv_types = None
    #     return self.sv_types
    #
    # def get_merge_count(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.vcf_merge_count = None
    #     if self.extension.endswith(".vcf"):
    #         self.vcf_merge_count = payload.vcf_merge_count
    #     return self.vcf_merge_count
    #
    # def get_variant_lines(self, payload):
    #     if self.extension.endswith(".gvf"):
    #         self.variant_lines = payload.gvf_line_count
    #     if self.extension.endswith(".vcf"):
    #         self.variant_lines = payload.vcf_line_count
    #     return self.variant_lines
    #
    #
    # def get_sample_number_missing(self, samples, vcf_header_fields):
    #     if self.extension.endswith(".gvf"):
    #         self.sample_number_missing = None
    #     if self.extension.endswith(".vcf"):
    #         samples_in_vcf_header = len(vcf_header_fields.split("\t")) - 8  # eight mandatory VCF header fields
    #         self.sample_number_missing = abs(samples_in_vcf_header - len(samples))
    #         #TODO: number of times a sample is missing in the merged vcf line, count per sample
    #     return self.sample_number_missing
    #
    # def get_sample_number(self, samples, vcf_header_fields):
    #     # if GVF, get the output from convert_gvf_pragmas_for_vcf_header
    #     if self.extension.endswith(".gvf"):
    #         self.sample_number = len(samples)
    #     # if VCF, get the header fields from generate_vcf_header_line
    #     if self.extension.endswith(".vcf"):
    #         header_fields = vcf_header_fields.split("\t")
    #         if "FORMAT" not in header_fields:
    #             self.sample_number = len(header_fields) - 8 # eight mandatory VCF header fields
    #         else:
    #             self.sample_number = len(header_fields) - 9 # 9 = eight mandatory VCF header fields + FORMAT
    #     return self.sample_number

    @staticmethod
    def get_file_md5(path_to_file):
        hash_md5 = hashlib.md5()
        with open(path_to_file, "rb") as f:
            # Read file in 4KB chunks until you get to an empty byte string
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        md5 = hash_md5.hexdigest()
        return md5

    @staticmethod
    def find_version(search_term, header_list):
        for line in header_list:
            if search_term in line:
                return line
        return None

    # def get_file_version(self, file_header):
    #     if self.extension.endswith(".gvf"):
    #         # this comes from the pragma that starts with: ##gvf-version 1.06
    #         self.version = self.find_version("##gvf-version", file_header)
    #     elif self.extension.endswith(".vcf"):
    #         # this comes from the VCF header that starts with: ##fileformat=VCFv4.4
    #         self.version = self.find_version("##fileformat", file_header)
    #     else:
    #         raise ValueError(f"Unsupported file extension: {self.extension}. Expected values: .gvf or .vcf")
    #     return self.version

    def __str__(self):
        summary_string = (f"======GVF======\n"
                          f"File name = {self.gvf_file_name}\n"
                          f"File path = {self.gvf_file_path}\n"
                          f"File extension = {self.gvf_extension}\n"
                          f"File size = {self.gvf_file_size} bytes\n"
                          f"Last mod Timestamp = {self.gvf_last_modified_timestamp}\n"
                          f"Version = {self.gvf_version}\n"
                          f"md5 string = {self.gvf_md5}\n"
                          f"gvf chrom = {self.gvf_chromosome_count}\n"
                          f"vcf chrom = {self.vcf_chromosome_count}\n"
                          f"gvf_feature_line_count = {self.gvf_feature_line_count}\n"
                          f"vcf_data_line_count = {self.vcf_data_line_count}\n"
                          f"vcf_number_of_merges = {self.vcf_number_of_merges}\n"
                          f"vcf_alt_allele_count = {self.vcf_alt_alleles_count}\n"
                          f"vcf_info_counter = {self.vcf_info_counter}\n"
                          f"=====\n"
                          )
        return summary_string
