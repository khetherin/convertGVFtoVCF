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

    def test_aggregation(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        expected_submitterDetails = [
            {'lastName': 'Doe', 'firstName': 'Jane', 'email': 'jane@example.com', 'laboratory': 'Lab B', 'centre': 'Institute B'},
            {'lastName': 'Smith', 'firstName': 'John', 'email': 'john@example.com', 'laboratory': 'Lab A', 'centre': 'Institute A'}
        ]
        expected_project = {
            'title': 'Project Alpha', 'description': 'Study A', 'taxId': 9606, 'centre': 'Institute A'
        }
        expected_analysis = [
            {'analysisTitle': 'Analysis B', 'analysisAlias': 'b2', 'description': 'Desc B', 'experimentType': 'Exome sequencing', 'referenceGenome': 'GCA_000001635.9'},
            {'analysisTitle': 'Analysis A', 'analysisAlias': 'a1', 'description': 'Desc A', 'experimentType': 'Whole genome sequencing', 'referenceGenome': 'GCA_000001405.15'}
        ]
        expected_sample = [
            {'analysisAlias': ['b2'], 'sampleInVCF': 'Sample_02', 'bioSampleAccession': 'SAMN00000002'},
            {'analysisAlias': ['a1'], 'sampleInVCF': 'Sample_01', 'bioSampleAccession': 'SAMN00000001'}
        ]
        expected_files = [
            {'analysisAlias': 'b2', 'fileName': 'data2.vcf.gz'},
            {'analysisAlias': 'a1', 'fileName': 'data1.vcf.gz'}
        ]
        metadata_from_file_1 = {
            "submitterDetails": [expected_submitterDetails[0]],
            "project": expected_project,
            "analysis": [expected_analysis[0]],
            "sample": [expected_sample[0]],
            "files": [expected_files[0]]
        }
        metadata_from_file_2 = {
            "submitterDetails": [expected_submitterDetails[1]],
            "project": expected_project,
            "analysis": [expected_analysis[1]],
            "sample": [expected_sample[1]],
            "files": [expected_files[1]]
        }
        json_objects_list = [metadata_from_file_1, metadata_from_file_2]
        aggregated_result = {
            "submitterDetails": coordinator.aggregate_submitter_details(json_objects_list),
            "project": coordinator.aggregate_project(json_objects_list),
            "analysis": coordinator.aggregate_analysis(json_objects_list),
            "sample": coordinator.aggregate_sample(json_objects_list),
            "files": coordinator.aggregate_files(json_objects_list)
        }
        self.assertEqual(aggregated_result["submitterDetails"], expected_submitterDetails)
        self.assertEqual(aggregated_result["project"], expected_project)
        self.assertEqual(aggregated_result["analysis"], expected_analysis)
        self.assertEqual(aggregated_result["sample"], expected_sample)
        self.assertEqual(aggregated_result["files"], expected_files)

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

    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.eva_update_metadata_with_vcf')
    @patch('convert_gvf_to_vcf.gvf_metadata_coordinator.convert')
    def test_convert_individual_gvf(self, mock_convert, mock_update):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        coordinator.base_output_dir = "tests/output"
        coordinator.project_paths = self.paths

        study_accession = "estd1"
        vcf_file = "tests/output/submission/estd1_Redon_et_al_2006/estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.vcf"
        input_gvf = "tests/data_dir/estd1_Redon_et_al_2006/gvf/estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf"

        coordinator.convert_individual_gvf(
            assembly_path="tests/input/human.fa",
            eva_retriever="mock_retriever",
            individual_gvf=input_gvf,
            json_eva="eva.json"
        )

        mock_convert.assert_called_once_with(
            gvf_input=input_gvf,
            vcf_output=vcf_file,
            assembly="tests/input/human.fa",
            paths=self.paths
        )
        mock_update.assert_called_once_with(
            eva_retriever="mock_retriever",
            json_eva="eva.json",
            vcf_output=vcf_file,
            study_accession=study_accession
        )

    def test_reconfigure_metadata(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        mock_eva_retriever = MagicMock()

        json_eva = "eva_metadata.json"
        study_accession = "estd1"
        study_name = "estd1_Redon_et_al_2006"
        date = "2014-04-01"
        submitted_files = [f"{study_name}.{date}.GRCh37.Submitted.gvf", f"{study_name}.{date}.GRCh38.Submitted.gvf",]
        remapped_files = [f"{study_name}.{date}.GRCh37.Remapped.gvf", f"{study_name}.{date}.GRCh38.Remapped.gvf"]

        coordinator._determine_same_and_reconfigure_json = MagicMock()
        coordinator.reconfigure_metadata(
            eva_retriever=mock_eva_retriever,
            json_eva=json_eva,
            remapped_files=remapped_files,
            study_accession=study_accession,
            submitted_files=submitted_files
        )

        self.assertEqual(coordinator._determine_same_and_reconfigure_json.call_count, 2)

        coordinator._determine_same_and_reconfigure_json.assert_any_call(
            study_accession, submitted_files, json_eva, mock_eva_retriever
        )
        coordinator._determine_same_and_reconfigure_json.assert_any_call(
            study_accession, remapped_files, json_eva, mock_eva_retriever
        )


    def test_process_single_assembly(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date = "2014-04-01"
        files_in_assembly = [f"{study_name}.{date}.GRCh37.Remapped.gvf",
                     f"{study_name}.{date}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date}.GRCh37.Submitted.gvf"
        ]
        assembly_name = "GRCh38"
        study_accession = "estd1"
        gvf_name_groups = {
            files_in_assembly[0]: (study_name, date, assembly_name)
        }

        coordinator.set_up_inputs_and_outputs = MagicMock(return_value=(
            "/path/to/assembly", "/path/to/report", "dgva.json", "eva.json", "out.vcf"
        ))

        mock_eva_retriever = MagicMock()
        mock_dgva_retriever = MagicMock()
        coordinator.retrieve_metadata = MagicMock(return_value=(
            mock_eva_retriever, mock_dgva_retriever
        ))

        coordinator.check_files = MagicMock(return_value=(
            ["remapped.gvf"], ["submitted.gvf"]
        ))

        result = coordinator._process_single_assembly(
            assembly_name, files_in_assembly, gvf_name_groups, study_accession
        )

        assembly_path, eva_retriever, json_eva, remapped_files, submitted_files = result

        self.assertEqual(assembly_path, "/path/to/assembly")
        self.assertEqual(eva_retriever, mock_eva_retriever)
        self.assertEqual(json_eva, "eva.json")
        self.assertEqual(remapped_files, ["remapped.gvf"])
        self.assertEqual(submitted_files, ["submitted.gvf"])

        coordinator.set_up_inputs_and_outputs.assert_called_once_with(study_name,date, "GRCh38")
        coordinator.retrieve_metadata.assert_called_once_with(
            "eva.json", "dgva.json", study_accession, "/path/to/assembly", "/path/to/report"
        )
        coordinator.check_files.assert_called_once_with(files_in_assembly)

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
                     f"{study_name}.{date}.GRCh37.Submitted.gvf"
        ]
        study_1, date_1, assembly_1 = coordinator.parse_gvf_filename(gvf_files[0])
        study_2, date_2, assembly_2 = coordinator.parse_gvf_filename(gvf_files[1])
        study_3, date_3, assembly_3 = coordinator.parse_gvf_filename(gvf_files[2])

        self.assertEqual(study_1, study_name)
        self.assertEqual(date_1, date)
        self.assertEqual(assembly_1, "GRCh37")

        self.assertEqual(study_2, study_name)
        self.assertEqual(date_2, date)
        self.assertEqual(assembly_2, "GRCh38")

        self.assertEqual(study_3, study_name)
        self.assertEqual(date_3, date)
        self.assertEqual(assembly_3, "GRCh37")

    def test_determine_same_and_reconfigure_json(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_accession = "estd1"
        study_name = "estd1_Redon_et_al_2006"
        date_1 = "2014-04-01"
        date_2 = "2015-04-01"
        gvf_files = [f"{study_name}.{date_1}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date_2}.GRCh38.Remapped.gvf"
        ]
        json_eva = "eva.json"
        mock_eva_retriever = MagicMock()
        mock_eva_retriever._fetch_analysis_reference_genome.return_value = "GRCh38"
        mock_eva_retriever._fetch_analysis_analysis_type.return_value = ["Read depth and paired-end mapping"]
        mock_eva_retriever._fetch_analysis_method_type.return_value = ["Sequencing"]
        mock_eva_retriever._determine_analysis_experiment_type.return_value = ["Whole genome sequencing"]

        coordinator._reconfigure_json_multi_analysis = MagicMock()

        coordinator._determine_same_and_reconfigure_json(
            study_accession=study_accession,
            submitted_or_remapped_files=gvf_files,
            json_eva=json_eva,
            eva_retriever=mock_eva_retriever
        )
        self.assertEqual(mock_eva_retriever._fetch_analysis_reference_genome.call_count, 2)
        coordinator._reconfigure_json_multi_analysis.assert_called_once_with(
            json_eva, gvf_files
        )

    def test_reconfigure_json_multi_analysis(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date_1 = "2014-04-01"
        date_2 = "2015-04-01"
        gvf_files = [f"{study_name}.{date_1}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date_2}.GRCh38.Remapped.gvf"
        ]

        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.json') as temp_file:
            json.dump({
                "submitterDetails": [],
                "project": {},
                "analysis": [{}],
                "sample": [{}],
                "files": [{}]
            }, temp_file)
            temp_json_path = temp_file.name
        coordinator._update_analysis_and_file_blocks = MagicMock(return_value=(
                ["updated_analysis"], ["updated_files"], ["new_alias"]
            ))
        coordinator._update_sample_block = MagicMock()


        parsed_output = coordinator._reconfigure_json_multi_analysis(temp_json_path, gvf_files)

        coordinator._update_analysis_and_file_blocks.assert_called_once_with([{}], gvf_files, [{}])
        coordinator._update_sample_block.assert_called_once_with(["new_alias"], [{}])

        self.assertEqual(parsed_output["analysis"], ["updated_analysis"])
        self.assertEqual(parsed_output["files"], ["updated_files"])
        self.assertEqual(parsed_output["sample"], [{}])

    def test_update_analysis_and_file_blocks(self):
        coordinator = GvfMetadataCoordinator(MagicMock(), MagicMock(), MagicMock())
        study_name = "estd1_Redon_et_al_2006"
        date_1 = "2014-04-01"
        date_2 = "2015-04-01"
        gvf_files = [f"{study_name}.{date_1}.GRCh38.Remapped.gvf",
                     f"{study_name}.{date_2}.GRCh38.Remapped.gvf"
        ]
        analysis_list = [{"analysisTitle": "ProjectOne", "analysisAlias": "study1_alias"}]
        files_list = [{"fileType": "vcf", "fileName": "old.vcf", "analysisAlias": "study1_alias"}]
        result = coordinator._update_analysis_and_file_blocks(
            analysis_list, gvf_files, files_list
        )
        multiple_analyses, multiple_files, new_analysis_aliases = result
        expected_aliases = ["study1_alias_file_1", "study1_alias_file_2"]
        self.assertEqual(new_analysis_aliases, expected_aliases)

        self.assertEqual(len(multiple_analyses), 2)
        self.assertEqual(multiple_analyses[0]["analysisAlias"], "study1_alias_file_1")
        self.assertEqual(multiple_analyses[0]["analysisTitle"], "ProjectOne")
        self.assertEqual(multiple_analyses[1]["analysisAlias"], "study1_alias_file_2")

        self.assertEqual(len(multiple_files), 2)
        self.assertEqual(multiple_files[0]["analysisAlias"], "study1_alias_file_1")
        self.assertEqual(multiple_files[0]["fileName"], f"{study_name}.{date_1}.GRCh38.Remapped.gvf")
        self.assertEqual(multiple_files[0]["fileType"], "vcf")

        self.assertEqual(multiple_files[1]["analysisAlias"], "study1_alias_file_2")
        self.assertEqual(multiple_files[1]["fileName"], f"{study_name}.{date_2}.GRCh38.Remapped.gvf")

if __name__ == '__main__':
    unittest.main()
