import os
import json
import tempfile
from unittest import TestCase
from unittest.mock import patch, MagicMock
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

        mock_retriever = MockRetriever.return_value
        mock_retriever._get_file_name.return_value = "test.vcf"
        mock_retriever._get_file_size.return_value = 1024
        mock_retriever._get_file_md5.return_value = "abc123md5"

        with tempfile.TemporaryDirectory() as tmpdir:
            json_output = os.path.join(tmpdir, "output.json")
            vcf_output = os.path.join(tmpdir, "test.vcf")

            initial_data = {
                "analysis": [{"analysisAlias": "alias_1"}],
                "files": []
            }
            with open(json_output, "w") as f:
                json.dump(initial_data, f)

            eva_add_file_metadata(mock_retriever, json_output, vcf_output, "estd1")

            expected_backup_path = os.path.join(tmpdir, "output_preconverted.json")
            self.assertTrue(os.path.exists(expected_backup_path))

            with open(json_output, "r") as f:
                updated_metadata = json.load(f)

            target_file = updated_metadata["files"][0]
            self.assertEqual(target_file["fileName"], "test.vcf")
            self.assertEqual(target_file["fileSize"], 1024)
            self.assertEqual(target_file["md5"], "abc123md5")
            self.assertEqual(target_file["analysisAlias"], "alias_1")
