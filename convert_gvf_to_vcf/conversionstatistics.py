"""
The purpose of this file is to calculate the (conversion) statistics of the GVF and VCF file.
"""
import os
from dataclasses import dataclass
from datetime import datetime
import hashlib
from typing import Optional, Union


@dataclass
class StatisticsPayload:
    """ Data transfer object for statistics input"""
    samples: list[str]
    vcf_header_fields: Optional[Union[str, None]]
    gvf_line_count: Optional[Union[int, None]]
    vcf_line_count: Optional[Union[int, None]]
    gvf_chromosome_count: Optional[Union[dict, None]]
    vcf_chromosome_count: Optional[Union[dict, None]]

class FileStatistics:
    """
    The responsibility of this class is to determine the statistics of a GVF or VCF file.
    """
    def __init__(self, file_path, payload: StatisticsPayload):
        self.file_path = file_path
        # file name
        self.file_name = os.path.basename(self.file_path)
        self.extension = os.path.splitext(self.file_path)[1]
        # file metadata
        self.metadata_stats = os.stat(self.file_path)
        self.file_size = self.metadata_stats.st_size
        last_modified_timestamp_raw = self.metadata_stats.st_mtime
        self.last_modified_timestamp = datetime.fromtimestamp(last_modified_timestamp_raw)
        self.version = None
        self.md5 = self.get_file_md5()
        # biological samples
        self.sample_number = self.get_sample_number(payload.samples, payload.vcf_header_fields)
        self.sample_number_missing = self.get_sample_number_missing(payload.samples, payload.vcf_header_fields)
        # chromsomes
        self.gvf_chromosome_naming_convention = payload.gvf_chromosome_count
        self.vcf_chromosome_naming_convention = payload.vcf_chromosome_count
        # variants
        # TODO: number of lines for GVF/VCF
        # self.variant_lines = self.get_variant_lines()
        # TODO: number of additional ALT alleles
        # TODO: number of ALT alleles
        # TODO: number of ALT alleles that could not be found
        # TODO: number of SV types
        # TODO: number of precise variants
        # TODO: number of imprecise variants
        # attribute mapping
        # TODO: observed GVF attributes (name and count)
        # TODO: INFO FORMAT tags (name and count)
        # TODO: dropped GVF attributes (assistingconverter.py:convert_gvf_attributes_to_vcf_values catching the attributes for a later date)
        # TODO: number of SO terms


    def get_variant_lines(self):
        # TODO: count GVF lines, VCF lines and merge events
        pass

    def get_sample_number_missing(self, samples, vcf_header_fields):
        if self.extension.endswith(".gvf"):
            self.sample_number_missing = None
        if self.extension.endswith(".vcf"):
            samples_in_vcf_header = len(vcf_header_fields.split("\t")) - 8  # eight mandatory VCF header fields
            self.sample_number_missing = abs(samples_in_vcf_header - len(samples))
            #TODO: number of times a sample is missing in the merged vcf line, count per sample
        return self.sample_number_missing

    def get_sample_number(self, samples, vcf_header_fields):
        # if GVF, get the output from convert_gvf_pragmas_for_vcf_header
        if self.extension.endswith(".gvf"):
            self.sample_number = len(samples)
        # if VCF, get the header fields from generate_vcf_header_line
        if self.extension.endswith(".vcf"):
            header_fields = vcf_header_fields.split("\t")
            if "FORMAT" not in header_fields:
                self.sample_number = len(header_fields) - 8 # eight mandatory VCF header fields
            else:
                self.sample_number = len(header_fields) - 9 # 9 = eight mandatory VCF header fields + FORMAT
        return self.sample_number

    def get_file_md5(self):
        hash_md5 = hashlib.md5()
        with open(self.file_path, "rb") as f:
            # Read file in 4KB chunks until you get to an empty byte string
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        self.md5 = hash_md5.hexdigest()
        return self.md5

    def find_version(self, search_term, header_list):
        for line in header_list:
            if search_term in line:
                return line
        return None

    def get_file_version(self, file_header):
        if self.extension.endswith(".gvf"):
            # this comes from the pragma that starts with: ##gvf-version 1.06
            self.version = self.find_version("##gvf-version", file_header)
        elif self.extension.endswith(".vcf"):
            # this comes from the VCF header that starts with: ##fileformat=VCFv4.4
            self.version = self.find_version("##fileformat", file_header)
        else:
            raise ValueError(f"Unsupported file extension: {self.extension}. Expected values: .gvf or .vcf")
        return self.version

    def __str__(self):
        summary_string = (f"File name = {self.file_name}\n"
                          f"File path = {self.file_path}\n"
                          f"File extension = {self.extension}\n"
                          f"File size = {self.file_size} bytes\n"
                          f"Last mod Timestamp = {self.last_modified_timestamp}\n"
                          f"Version = {self.version}\n"
                          f"md5 string = {self.md5}\n"
                          f"=====\n"
                          f""
                          f"sample_numer = {self.sample_number}\n"
                          f"GVF chromosomes = {self.gvf_chromosome_naming_convention}\n"
                          f"VCF chromosomes = {self.vcf_chromosome_naming_convention}\n"
                          )
        return summary_string
