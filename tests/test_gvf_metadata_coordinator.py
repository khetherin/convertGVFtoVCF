import builtins
import json
import os.path
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open

from convert_gvf_to_vcf.gvf_metadata_coordinator import GvfMetadataCoordinator
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestGvfMetadataCoordinator(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.config = self.paths.full_config_path
        self.test_dir = self.paths.test_dir
        self.etc_folder = self.paths.etc_dir
        self.output_dir = os.path.join(self.test_dir, "output", "gvf_metadata_coord")
        os.makedirs(self.output_dir, exist_ok=True)
        self.paths = ProjectPaths()
        self.input_dir = Path(os.path.join(self.paths.test_dir, "input", "json_aggregator"))
        self.output_file = Path(os.path.join(self.paths.test_dir, "output", "json_aggregator", "eva_submission.json"))

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.logger')
    def test_process_studies(self, mock_logger):
        no_gvf_files = []
        single_gvf_file = ["path/to/file1.gvf"]
        multiple_gvf_files = ["path/to/file1.gvf", "path/to/file2.gvf"]

        # no gvf
        study_and_gvf_files_input_data = {"estd1": no_gvf_files}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator.base_output_dir = self.output_dir
        coordinator._process_no_gvf_files = MagicMock()
        coordinator.process_studies()
        coordinator._process_no_gvf_files.assert_called_once_with([], "estd1")
        # single gvf
        study_and_gvf_files_input_data = {"estd1": single_gvf_file}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator.base_output_dir = self.output_dir
        coordinator.parse_gvf_filename = MagicMock(return_value=("estd1_Redon_et_al_2006", None, None))
        coordinator._process_gvf_files = MagicMock()
        coordinator.process_studies()
        expected_json_1 = os.path.join(self.output_dir, "submission", "estd1_Redon_et_al_2006",
                                       "eva_submission_estd1.json")
        coordinator._process_gvf_files.assert_called_once_with(single_gvf_file, "estd1", expected_json_1)

        # multiple
        study_and_gvf_files_input_data = {"estd1": multiple_gvf_files}
        coordinator = GvfMetadataCoordinator(study_and_gvf_files_input_data, self.output_dir, self.config)
        coordinator.base_output_dir = self.output_dir
        coordinator.parse_gvf_filename = MagicMock(return_value=("estd1_Redon_et_al_2006", None, None))
        coordinator._process_gvf_files = MagicMock()
        coordinator.process_studies()
        expected_json_2 = os.path.join(self.output_dir, "submission", "estd1_Redon_et_al_2006",
                                       "eva_submission_estd1.json")
        coordinator._process_gvf_files.assert_called_once_with(multiple_gvf_files, "estd1", expected_json_2)

    def test_process_gvf_files(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())

        study_name = "estd1_Redon_et_al_2006"
        date = "2014-04-01"
        assembly_name = "GRCh38"
        study_accession = "estd1"
        master_json = "tests/output/master_metadata.json"

        files_in_assembly = [f"{study_name}.{date}.{assembly_name}.Remapped.gvf"]

        mock_assembly_to_paths = {assembly_name: files_in_assembly}
        mock_gvf_name_groups = {files_in_assembly[0]: (study_name, date, assembly_name)}

        coordinator._group_files_by_assembly = MagicMock(
            return_value=(mock_assembly_to_paths, mock_gvf_name_groups))

        mock_initial_master = {
            "submitterDetails": None, "project": None, "analysis": [], "sample": [], "files": []
        }
        coordinator._initialize_master_metadata = MagicMock(return_value=mock_initial_master)
        coordinator.set_up_inputs_and_outputs = MagicMock(return_value=(
            "/path/to/assembly", "/path/to/report", "dgva.json", "eva.json", "out.vcf"
        ))
        mock_eva_metadata = {
            "submitterDetails": [{"firstName": "John", "lastName": "Doe", "email": "john@doe.com"}],
            "project": {"title": "Test Variant Project"}
        }
        coordinator._stage_assembly_metadata = MagicMock(return_value=mock_eva_metadata)
        coordinator.determine_submitted_and_remapped_files = MagicMock(return_value=(
            ["remapped.gvf"], ["submitted.gvf"]
        ))
        coordinator.convert_individual_gvf = MagicMock(return_value="/path/to/submission/out.vcf")
        coordinator._add_multiple_assemblies_into_metadata = MagicMock()
        coordinator._save_master_metadata = MagicMock()
        coordinator._process_gvf_files(
            gvf_files=files_in_assembly,
            study_accession=study_accession,
            master_json=master_json
        )

        self.assertEqual(coordinator.assembly_path, "/path/to/assembly")
        self.assertEqual(coordinator.assembly_report_path, "/path/to/report")
        self.assertEqual(coordinator.json_dgva_path, "dgva.json")
        self.assertEqual(coordinator.json_eva_path, "eva.json")
        self.assertEqual(coordinator.vcf_output_path, "out.vcf")

        coordinator._group_files_by_assembly.assert_called_once_with(files_in_assembly)
        coordinator._initialize_master_metadata.assert_called_once()
        coordinator.set_up_inputs_and_outputs.assert_called_once_with(study_name, date, assembly_name)
        coordinator._stage_assembly_metadata.assert_called_once_with(assembly_name, study_accession)
        coordinator.determine_submitted_and_remapped_files.assert_called_once_with(files_in_assembly)
        coordinator.convert_individual_gvf.assert_called_once_with(files_in_assembly[0])
        coordinator._add_multiple_assemblies_into_metadata.assert_called_once_with(
            json_eva="eva.json",
            gvf_files=["remapped.gvf"],  # Uses remapped_files if populated
            master_metadata=mock_initial_master,
            vcf_path="/path/to/submission/out.vcf"
        )
        coordinator._save_master_metadata.assert_called_once_with(mock_initial_master, master_json)

    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.eva_update_metadata_with_vcf')
    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.convert')
    def test_convert_individual_gvf(self, mock_convert, mock_update):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        coordinator.base_output_dir = "tests/output"
        coordinator.project_paths = self.paths

        assembly = "tests/input/human.fa"
        coordinator.assembly_path = assembly
        coordinator.json_eva_path = "eva.json"
        coordinator.eva_retriever = "mock_retriever"

        study_accession = "estd1"
        vcf_file = "tests/output/submission/estd1_Redon_et_al_2006/estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.vcf"
        input_gvf = "tests/data_dir/estd1_Redon_et_al_2006/gvf/estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf"

        output_result = coordinator.convert_individual_gvf(individual_gvf=input_gvf)

        mock_convert.assert_called_once_with(
            gvf_input=input_gvf,
            vcf_output=vcf_file,
            assembly=assembly,
            paths=self.paths
        )

        mock_update.assert_called_once_with(
            retriever="mock_retriever",
            json_output="eva.json",
            vcf_output=vcf_file,
            study_accession=study_accession
        )

        self.assertEqual(output_result, vcf_file)

    def test_suffix_assembly_metadata(self):
        mock_metadata = {
            "analysis": [{"analysisAlias": "estd1_analysis", "analysisTitle": "estd1"}],
            "files": [{"analysisAlias": "estd1_analysis", "fileName": "test.vcf"}],
            "sample": [{"analysisAlias": ["estd1_analysis"], "sampleInVCF": "NA18570"}]
        }

        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        coordinator._suffix_assembly_metadata(mock_metadata, "GRCh38")


        self.assertEqual(mock_metadata["analysis"][0]["analysisAlias"], "estd1_analysis_GRCh38")
        self.assertEqual(mock_metadata["analysis"][0]["analysisTitle"], "estd1_GRCh38")
        self.assertEqual(mock_metadata["files"][0]["analysisAlias"], "estd1_analysis_GRCh38")
        self.assertEqual(mock_metadata["sample"][0]["analysisAlias"], ["estd1_analysis_GRCh38"])

    def test_cleanup_master_metadata(self):
        mock_master_metadata = {
            "sample": [
                {"sampleInVCF": "NA18570", "analysisAlias": ["estd1_GRCh38"]},
                {"sampleInVCF": "NA18570", "analysisAlias": ["estd1_NCBI35"]}
            ],
            "files": [
                {"analysisAlias": "estd1_NCBI35", "fileName": "submitted.vcf"},
                {"analysisAlias": "estd1_GRCh38", "fileName": "remapped.vcf"}
            ]
        }

        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())

        coordinator._cleanup_master_metadata(mock_master_metadata)
        sample = mock_master_metadata["sample"][0]["sampleInVCF"]
        analysisAlias = mock_master_metadata["sample"][0]["analysisAlias"]

        self.assertEqual(len(mock_master_metadata["sample"]), 1) # EXPECT ONE SAMPLE
        self.assertEqual(sample, "NA18570")

        self.assertEqual(len(analysisAlias), 2)
        self.assertIn("estd1_GRCh38", analysisAlias)
        self.assertIn("estd1_NCBI35", analysisAlias)

        self.assertEqual(len(mock_master_metadata["files"]), 2)
        self.assertEqual(mock_master_metadata["files"][0]["fileName"], "submitted.vcf")
        self.assertEqual(mock_master_metadata["files"][1]["fileName"], "remapped.vcf")

    def test_stage_assembly_metadata(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())

        assembly_name = "GRCh38"
        study_accession = "estd1"

        coordinator.assembly_path = "/path/to/assembly"
        coordinator.assembly_report_path = "/path/to/report"
        coordinator.json_dgva_path = "dgva.json"
        coordinator.json_eva_path = "eva.json"

        mock_eva_retriever = MagicMock()
        mock_dgva_retriever = MagicMock()
        coordinator.retrieve_metadata = MagicMock(return_value=(
            mock_eva_retriever, mock_dgva_retriever
        ))

        mock_initial_metadata = {"submitterDetails": "initial_data"}
        coordinator._load_json_file = MagicMock(return_value=mock_initial_metadata)
        coordinator._suffix_assembly_metadata = MagicMock()

        with patch.object(builtins, 'open', mock_open()):
            result_eva_metadata = coordinator._stage_assembly_metadata(assembly_name, study_accession)

        self.assertEqual(result_eva_metadata, mock_initial_metadata)
        coordinator.retrieve_metadata.assert_called_once_with(
            "eva.json",
            "dgva.json",
            study_accession,
            "/path/to/assembly",
            "/path/to/report"
        )
        coordinator._load_json_file.assert_called_once_with("eva.json")
        coordinator._suffix_assembly_metadata.assert_called_once_with(mock_initial_metadata, assembly_name)

    def test_group_files_by_assembly(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date = "2014-04-01"
        gvf_files = [f"{study_name}.{date}.GRCh37.Remapped.gvf",
                     f"{study_name}.{date}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date}.GRCh37.Submitted.gvf"
        ]

        coordinator.parse_gvf_filename = MagicMock(side_effect=[
            (study_name, date, "GRCh37"),
            (study_name, date, "GRCh38"),
            (study_name, date, "GRCh37")
        ])
        assembly_groups, gvf_name_groups = coordinator._group_files_by_assembly(gvf_files)

        self.assertEqual(assembly_groups["GRCh37"], [f"{study_name}.{date}.GRCh37.Remapped.gvf", f"{study_name}.{date}.GRCh37.Submitted.gvf"])
        self.assertEqual(assembly_groups["GRCh38"], [f"{study_name}.{date}.GRCh38.Remapped.gvf"])

        self.assertEqual(gvf_name_groups[f"{study_name}.{date}.GRCh37.Remapped.gvf"], (study_name, date, "GRCh37"))
        self.assertEqual(gvf_name_groups[f"{study_name}.{date}.GRCh38.Remapped.gvf"], (study_name, date, "GRCh38"))

    def test_parse_gvf_filename(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date = "2014-04-01"
        gvf_files = [f"{study_name}.{date}.GRCh37.Remapped.gvf",
                     f"{study_name}.{date}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date}.GRCh37.Submitted.gvf",
                     f"{study_name}.{date}.GRCh37.Submitted.gvf.sorted.gvf",
                     f"{study_name}.{date}.GRCh37.p13.Submitted.gvf",
                     f"{study_name}.{date}.GRCh37.p13.Submitted.gvf.sorted.gvf",
                     f"{study_name}.{date}.GRCh37.ABCDEFG.txt", # invalid

        ]
        study_1, date_1, assembly_1 = coordinator.parse_gvf_filename(gvf_files[0])
        study_2, date_2, assembly_2 = coordinator.parse_gvf_filename(gvf_files[1])
        study_3, date_3, assembly_3 = coordinator.parse_gvf_filename(gvf_files[2])
        study_4, date_4, assembly_4 = coordinator.parse_gvf_filename(gvf_files[3])
        study_5, date_5, assembly_5 = coordinator.parse_gvf_filename(gvf_files[4])
        study_6, date_6, assembly_6 = coordinator.parse_gvf_filename(gvf_files[5])
        study_7, date_7, assembly_7 = coordinator.parse_gvf_filename(gvf_files[6])

        self.assertEqual(study_1, study_name)
        self.assertEqual(date_1, date)
        self.assertEqual(assembly_1, "GRCh37")

        self.assertEqual(study_2, study_name)
        self.assertEqual(date_2, date)
        self.assertEqual(assembly_2, "GRCh38")

        self.assertEqual(study_3, study_name)
        self.assertEqual(date_3, date)
        self.assertEqual(assembly_3, "GRCh37")

        self.assertEqual(study_4, study_name)
        self.assertEqual(date_4, date)
        self.assertEqual(assembly_4, "GRCh37")

        self.assertEqual(study_5, study_name)
        self.assertEqual(date_5, date)
        self.assertEqual(assembly_5, "GRCh37.p13")

        self.assertEqual(study_6, study_name)
        self.assertEqual(date_6, date)
        self.assertEqual(assembly_6, "GRCh37.p13")

        self.assertEqual(study_7, None)
        self.assertEqual(date_7, None)
        self.assertEqual(assembly_7, None)

    def test_add_multiple_assemblies_into_metadata(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        gvf_files = [
            f"{study_name}.2014-04-01.GRCh38.Remapped.gvf",
            f"{study_name}.2015-04-01.GRCh38.Remapped.gvf"
        ]

        metadata_to_add = {
            "submitterDetails": [{"firstName": "NameTest"}],
            "project": {"projectAccession": "PRJEB12345"},
            "analysis": [{"analysisAlias": "AliasTest"}],
            "sample": [{"analysisAlias": ["AliasTest"], "sampleInVCF": "Sample1"}],
            "files": [{"analysisAlias": "AliasTest", "fileName": "test.vcf"}]
        }

        master_metadata = {
            "submitterDetails": None, "project": None, "analysis": [], "sample": [], "files": []
        }
        vcf_path = "/a/path/to/my.vcf"


        with patch.object(GvfMetadataCoordinator, '_update_analysis_and_file_blocks',
                          return_value=["new_alias"]) as mock_analysis_block, \
                patch.object(GvfMetadataCoordinator, '_update_sample_block') as mock_sample_block:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_json_path = os.path.join(temp_dir, "eva_metadata.json")
                with open(temp_json_path, 'w', encoding='utf-8') as temp_out:
                    json.dump(metadata_to_add, temp_out)

                result_metadata = coordinator._add_multiple_assemblies_into_metadata(
                    json_eva=temp_json_path,
                    gvf_files=gvf_files,
                    master_metadata=master_metadata,
                    vcf_path=vcf_path
                )

                self.assertEqual(result_metadata, metadata_to_add)

            mock_analysis_block.assert_called_once_with(
                metadata_to_add=metadata_to_add,
                files=gvf_files,
                master_metadata=master_metadata,
                vcf_path=vcf_path
            )

            mock_sample_block.assert_called_once_with(
                new_analysis_aliases=["new_alias"],
                metadata_to_add=metadata_to_add,
                master_metadata=master_metadata
            )

    def test_update_analysis_and_file_blocks(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date_1 = "2014-04-01"
        date_2 = "2015-04-01"
        gvf_files = [f"{study_name}.{date_1}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date_2}.GRCh38.Remapped.gvf"
        ]
        metadata_to_add = {
            "submitterDetails":[{"firstName": "J", "lastName": "D", "email": "a@b.com", "laboratory": "L", "centre": "C"}],
            "project":{"projectAccession": "PRJEB12345"},
            "analysis": [{"analysisTitle": "ProjectOne", "analysisAlias": "study1_alias"}],
            "sample":[{"analysisAlias": ["study1_alias"], "sampleInVCF": "S1", "bioSampleAccession": "SAMEA123"}],
            "files": [{"fileType": "vcf", "fileName": "old.vcf", "analysisAlias": "study1_alias"}]
        }
        master_metadata = {
            "submitterDetails": None,
            "project": None,
            "analysis": [],
            "sample": [],
            "files": []
        }
        vcf_path="/a/test/output/path.vcf"
        new_analysis_aliases = coordinator._update_analysis_and_file_blocks(
            metadata_to_add=metadata_to_add,
            files=gvf_files,
            master_metadata=master_metadata,
            vcf_path=vcf_path
        )
        expected_aliases = ["study1_alias_file_1", "study1_alias_file_2"]


        self.assertEqual(new_analysis_aliases, expected_aliases)

        self.assertEqual(len(master_metadata["analysis"]), 2)
        self.assertEqual(master_metadata["analysis"][0]["analysisAlias"], "study1_alias_file_1")
        self.assertEqual(master_metadata["analysis"][0]["analysisTitle"], "ProjectOne")
        self.assertEqual(master_metadata["analysis"][1]["analysisAlias"], "study1_alias_file_2")
        self.assertEqual(master_metadata["analysis"][1]["analysisTitle"], "ProjectOne")

        self.assertEqual(len(master_metadata["files"]), 2)
        self.assertEqual(master_metadata["files"][0]["analysisAlias"], "study1_alias_file_1")
        self.assertEqual(master_metadata["files"][0]["fileName"], f"/a/test/output/path.vcf")
        self.assertEqual(master_metadata["files"][0]["fileType"], "vcf")

        self.assertEqual(master_metadata["files"][1]["analysisAlias"], "study1_alias_file_2")
        self.assertEqual(master_metadata["files"][1]["fileName"], f"/a/test/output/path.vcf")

if __name__ == '__main__':
    unittest.main()
