import os
import json
from unittest import TestCase
from unittest.mock import patch
from convert_gvf_to_vcf.gather_metadata import eva_add_file_metadata
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestGatherMetadata(TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.tests_folder = self.paths.test_dir
        self.config = os.path.join(self.tests_folder, "input", "test.config")
        self.study_accession = "estd22"
        self.assembly = os.path.join(self.tests_folder, "input", "zebrafish.fa")
        self.assembly_report = os.path.join(self.tests_folder, "input", "assembly_report.txt")

        self.json_file_eva = os.path.join(self.tests_folder, "output", "a.json")
        self.json_file_eva_preconverted = os.path.join(self.tests_folder, "output", "a_preconverted.json")
        self.vcf_output = os.path.join(self.tests_folder, "output", "a.vcf")
        self.expected_json_eva_output = os.path.join(self.tests_folder, "output", "a.json")
        self.json_file_dgva = os.path.join(self.tests_folder, "output", "a_dgva.json")


    @patch('convert_gvf_to_vcf.gather_metadata.EVAMetadataRetriever')
    def test_add_file_metadata(self, MockRetriever):
        self.json_file_eva = "tmp_eva.json"
        self.expected_json_eva_output = "tmp_expected.json"
        with open(self.json_file_eva, "w") as f:
            json.dump({"files": [{}]}, f)
        with open(self.expected_json_eva_output, "w") as f:
            json.dump({"files": [{"fileName": "a.vcf", "fileSize": 12345, "md5": "abcde12345"}]}, f)


        mock_retriever = MockRetriever.return_value
        mock_retriever._get_file_name.return_value = "a.vcf"
        mock_retriever._get_file_size.return_value = 12345
        mock_retriever._get_file_md5.return_value = "abcde12345"
        eva_add_file_metadata(mock_retriever, self.json_file_eva, self.vcf_output)

        with open(self.json_file_eva, 'r') as f:
            metadata = json.load(f)
        with open(self.expected_json_eva_output, 'r') as f_out:
            expected_metadata = json.load(f_out)
        self.assertEqual(metadata["files"][0]["fileName"], expected_metadata["files"][0]["fileName"])
        self.assertEqual(metadata["files"][0]["fileSize"], expected_metadata["files"][0]["fileSize"])
        self.assertEqual(metadata["files"][0]["md5"], expected_metadata["files"][0]["md5"])

        if os.path.exists(self.json_file_eva):
            os.remove(self.json_file_eva)
        if os.path.exists(self.expected_json_eva_output):
            os.remove(self.expected_json_eva_output)
