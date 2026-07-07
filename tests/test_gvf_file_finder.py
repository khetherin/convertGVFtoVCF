import unittest
import os
from unittest.mock import MagicMock, patch

from convert_gvf_to_vcf.gvf_file_finder import GvfFileFinder
from convert_gvf_to_vcf.project_paths import ProjectPaths


class TestGvfFileFinder(unittest.TestCase):
    def setUp(self):
        self.paths = ProjectPaths()
        self.top_dir = os.path.join(self.paths.test_dir, "data_dir")
        self.current_dir = os.path.join(self.top_dir, "estd1_Redon_et_al_2006", "gvf")
        self.study_accession = "estd1"

    def test_deduplicate_files(self):
        files = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.NCBI35.Submitted.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.p13.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh38.Remapped.gvf']
        expected_paths = []
        for file in files:
            expected_paths.append(os.path.join(self.current_dir, file))
        study_and_files = {self.study_accession: expected_paths}
        self._check_md5_match = MagicMock(return_value=False)
        os.path.getsize = MagicMock(return_value=100)
        result = GvfFileFinder.deduplicate_files(self, study_and_files)
        assert "estd1" in result
        assert result["estd1"] == expected_paths

    @patch('os.path.isdir', return_value=True)
    @patch('os.listdir')
    def test_find_study_dirs(self, mock_listdir, mock_isdir):
        file_finder = GvfFileFinder(search_dir=self.top_dir)

        mock_listdir.return_value = [
            'nstd1_Tuzun_et_al_2005',
            'estd1_Redon_et_al_2006',
            '.hidden_file',
            'wrong_prefix'
        ]

        result_all = file_finder._find_study_dirs(study_accession=None)
        self.assertEqual(result_all, ['estd1_Redon_et_al_2006', 'nstd1_Tuzun_et_al_2005'])

        result_specific = file_finder._find_study_dirs(study_accession="estd1")
        self.assertEqual(result_specific, ['estd1_Redon_et_al_2006'])

    @patch('os.path.isfile', return_value=True)
    @patch('os.listdir')
    def test_collect_gvf_files(self, mock_listdir, mock_isfile):
        file_finder = GvfFileFinder(search_dir=self.top_dir)
        gvf_dir = self.current_dir

        mock_listdir.return_value = [
            'file1.gvf',
            'file2.GVF',
            'file3.gvf.gz',
            'readme.txt',
            'image.png'
        ]
        result = file_finder._collect_gvf_files(gvf_dir)
        expected_paths = [
            os.path.join(gvf_dir, 'file1.gvf'),
            os.path.join(gvf_dir, 'file2.GVF')
        ]
        self.assertEqual(result, expected_paths)

    def test_scan(self):
        file_finder = GvfFileFinder(search_dir=self.top_dir)
        study_name = "estd1_Redon_et_al_2006"

        study_path = os.path.join(self.top_dir, study_name)
        gvf_path = os.path.join(study_path, "gvf")

        files = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.NCBI35.Submitted.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.p13.Remapped.gvf', 'estd1_Redon_et_al_2006.2014-04-01.GRCh38.Remapped.gvf']
        expected_paths = []
        for file in files:
            expected_paths.append(os.path.join(gvf_path, file))
        mock_return_dict = {self.study_accession: expected_paths} # assuming no duplicates paths
        file_finder._find_study_dirs = MagicMock(return_value=[study_name])
        file_finder._collect_gvf_files = MagicMock(return_value=files)
        file_finder.deduplicate_files = MagicMock(return_value=mock_return_dict)
        file_finder._process_gvf_directory = MagicMock(return_value =mock_return_dict)

        with patch('os.path.isdir', return_value=True):
            result = file_finder.scan(self.study_accession)

        self.assertEqual(result[self.study_accession], expected_paths)
        file_finder.deduplicate_files.assert_called_once()

    def test_get_md5(self):
        file_finder = GvfFileFinder(search_dir=self.top_dir)
        file_to_test = os.path.join(self.current_dir, 'estd1_Redon_et_al_2006.2014-04-01.GRCh37.Remapped.gvf')
        actual_md5 = file_finder.get_md5(file_to_test)
        expected_md5 = "0e1ae1a712ac90501f4a57ec90bff713"
        self.assertEqual(actual_md5, expected_md5)