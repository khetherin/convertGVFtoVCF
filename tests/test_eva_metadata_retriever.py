import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, Mock
from collections import namedtuple

from convert_gvf_to_vcf.metadata_retrievers.eva_metadata import EVAMetadataRetriever
from convert_gvf_to_vcf.project_paths import ProjectPaths

class TestEVAMetadataRetriever(TestCase):
    def setUp(self):
        input_folder = os.path.dirname(__file__)
        self.config = os.path.join(input_folder, "input", "test.config")
        self.vcf_output = os.path.join(input_folder, "output", "a.vcf")
        self.paths = ProjectPaths()
        self.full_config_path = self.paths.full_config_path

    def test_init(self):
        # testing the __init__ function has loaded the config correctly
        metadata_client = EVAMetadataRetriever(self.config)
        assert metadata_client._host == "mockhost"
        assert metadata_client._port == "1111"
        assert metadata_client._username == "mockuser"
        assert metadata_client._password == "mockpassword"
        assert metadata_client._service_name == "mockservicename"
        assert metadata_client._connection == None
        assert metadata_client._max_retries == 3

    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_connection_new(self, mock_connection):
        mock_connection_object = Mock()
        mock_connection.return_value = mock_connection_object
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client.connection
        self.assertEqual(result, mock_connection_object)
        mock_connection.assert_called_once()

    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_connection_healthy(self, mock_connection):
        # expected behaviour: if healthy, keep the connection
        mock_connection_object_existing = Mock()
        mock_connection.return_value = mock_connection_object_existing
        metadata_client = EVAMetadataRetriever(self.config)
        metadata_client._connection = mock_connection_object_existing
        result = metadata_client.connection
        self.assertEqual(result, mock_connection_object_existing)
        mock_connection.assert_not_called()

    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_connection_unhealthy(self, mock_connection):
        # expected behaviour: if unhealthy, close the connection then open a new connection
        mock_unhealthy = Mock()
        mock_unhealthy.is_healthy.return_value = False
        mock_new_connection = Mock()
        mock_connection.return_value = mock_new_connection

        metadata_client = EVAMetadataRetriever(self.config)
        # set connection to unhealthy
        metadata_client._connection = mock_unhealthy
        result = metadata_client.connection
        # did you close the unhealthy connection
        mock_unhealthy.close.assert_called_once()
        # did you get a new connection
        self.assertEqual(result, mock_new_connection)

    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_connection_max_retry_third_attempt_success(self, mock_connection):
        # expected behaviour: unsuccessful on first two attempts, successful on third attempt.
        mock_healthy = Mock()
        mock_healthy.is_healthy.return_value = True
        mock_connection.side_effect = [Exception, Exception, mock_healthy]

        metadata_client = EVAMetadataRetriever(self.config)
        # checking attribute is set to max 3 retries
        self.assertEqual(metadata_client._max_retries, 3)
        result = metadata_client.connection
        # checking your connection is healthy
        self.assertEqual(result, mock_healthy)
        self.assertEqual(mock_connection.call_count, 3)


    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_connection_max_retry_third_attempt_fail(self, mock_connection):
        # expected behaviour: it will allow 3 retries. it will not allow the fourth attempt
        mock_connection.side_effect = [Exception("attempt1"), Exception("attempt2"), Exception("attempt3"), Exception("attempt4")]
        metadata_client = EVAMetadataRetriever(self.config)
        with self.assertRaises(Exception) as unsuccessful:
            result = metadata_client.connection
        # this should only allow 3 attempts
        assert mock_connection.call_count == 3

    @patch("convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever.oracledb.connect")
    def test_load_from_db(self, mock_connect):
        mock_connection = MagicMock()
        mock_cursor = MagicMock()

        mock_connect.return_value = mock_connection
        mock_connection.cursor.return_value.__enter__.return_value = mock_cursor
        expected_data = [('row1',), ('row2',)]
        mock_cursor.fetchall.return_value = expected_data

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client.load_from_db("SQL_QUERY")
        self.assertEqual(result, expected_data)

    def test_create_json_file(self):
        pass

    def test__get_validated_value(self):
        pass

    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_addresses")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_centres")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_email_addresses")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_phone_numbers")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_first_names")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_submitter_details_all_last_names")
    def test_get_submitter_details(self,
                                   mock_fetch_all_last_names,
                                   mock_fetch_all_first_names,
                                   mock_fetch_all_phone_numbers,
                                   mock_fetch_all_email_addresses,
                                   mock_fetch_all_centres,
                                   mock_fetch_all_addresses):

        # first_name, last_name, telephone, email, laboratory, centre, address
        # create mock data from load_from_db
        mock_data_all_last_names = ["Smith", None, ""] # expected_input, None value, empty string
        mock_data_all_first_names = ["John", None, ""] # expected_input, None value, empty string
        mock_data_all_phone_numbers = ["0123456789", None, ""] # expected_input, None value, empty string
        mock_data_all_email_addresses = ["e1@mail.com", None, ""] # expected_input, None value, empty string
        mock_data_all_centres = ["centre1", None, ""] # expected_input, None value, empty string
        mock_data_all_addresses = ["123 Road Name, City, Country", None, ""] # expected_input, None value, empty string

        mock_fetch_all_last_names.return_value = mock_data_all_last_names
        mock_fetch_all_first_names.return_value = mock_data_all_first_names
        mock_fetch_all_phone_numbers.return_value = mock_data_all_phone_numbers
        mock_fetch_all_email_addresses.return_value = mock_data_all_email_addresses
        mock_fetch_all_centres.return_value = mock_data_all_centres
        mock_fetch_all_addresses.return_value = mock_data_all_addresses

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_submitter_details("STUDY123")

        # default value only, the None and empty string should only have empty strings for required fields
        expected = [{
                'firstName': 'John',
                'lastName': 'Smith',
                'telephone': '0123456789',
                'laboratory': '',
                'email': 'e1@mail.com',
                'centre': 'centre1',
                'address': '123 Road Name, City, Country'
            },
            {'lastName': '', 'firstName': '', 'email': '', 'laboratory': '', 'centre': ''},
            {'lastName': '', 'firstName': '', 'email': '', 'laboratory': '', 'centre': ''}
        ]
        self.assertEqual(expected, result)
        mock_fetch_all_last_names.assert_called_once_with("STUDY123")
        mock_fetch_all_first_names.assert_called_once_with("STUDY123")
        mock_fetch_all_phone_numbers.assert_called_once_with("STUDY123")
        mock_fetch_all_email_addresses.assert_called_once_with("STUDY123")
        mock_fetch_all_centres.assert_called_once_with("STUDY123")
        mock_fetch_all_addresses.assert_called_once_with("STUDY123")

    def test__get_project_pre_registered(self):
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_project_pre_registered("PRJEB123456789")
        expected = {'projectAccession': 'PRJEB123456789'}
        self.assertEqual(expected, result)

    # @patch("convert_gvf_to_vcf.metadataJSON.DGVaMetadataRetriever.validate_project")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_hold_date")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_project_links")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_project_parent_project")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_project_publications")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_centre")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_tax_id")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_project_description")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_project_title")
    def test__get_project_new(self, mock_title, mock_description, mock_tax_id, mock_centre, mock_publications, mock_parent, mock_links, mock_date):
        # all present
        mock_title.return_value = "title"
        mock_description.return_value = "description"
        mock_tax_id.return_value = int("9606")
        mock_centre.return_value = ["centre"]
        mock_publications.return_value = ['12345']
        mock_parent.return_value = "PRJNB5502"
        mock_links.return_value =  ["www.link.com", "www.another.com"]
        mock_date.return_value = str("2018-11-13")
        centre = mock_centre.return_value[0]
        publications = ['PubMed:12345']

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_project_new("STUDY123")
        expected = {"title": mock_title.return_value,
                    "description": mock_description.return_value,
                    "taxId": mock_tax_id.return_value,
                    "centre":centre,
                    "publications": publications,
                    "parentProject": mock_parent.return_value,
                    "links": mock_links.return_value,
                    "holdDate": mock_date.return_value
        }
        self.assertEqual(expected, result)
        # requried only
        mock_title.return_value = "title"
        mock_description.return_value = "description"
        mock_tax_id.return_value = int("9606")
        mock_centre.return_value = ["centre"]
        mock_publications.return_value = None
        mock_parent.return_value = None
        mock_links.return_value =  None
        mock_date.return_value = None
        centre = mock_centre.return_value[0]

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_project_new("STUDY123")
        expected = {"title": mock_title.return_value,
                    "description": mock_description.return_value,
                    "taxId": mock_tax_id.return_value,
                    "centre":centre
        }
        self.assertEqual(expected, result)

    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_run_accessions")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_pipeline_descriptions")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_software")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_platform")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_reference_genome")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._determine_analysis_experiment_type")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_method_type")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_analysis_type")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_description")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_alias")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_ids")
    def test__get_analysis(self,
                           mock_ids,
                           mock_alias,
                           mock_desc,
                           mock_analysis_type,
                           mock_method_type,
                           mock_determined_experiment_type,
                           mock_reference_genome,
                           mock_platform,
                           mock_software,
                           mock_pipeline_descriptions,
                           mock_run_accessions):
        # testing all values
        mock_ids.return_value = ["1","2", "3"]
        mock_alias.return_value = "MYanalysisALIAS_ALL"
        mock_desc.return_value = "mock_desc"
        mock_analysis_type.return_value = "Genotyping"
        mock_method_type.return_value = "Sequencing"
        mock_determined_experiment_type.return_value = "genotyping_by_sequencing"
        mock_reference_genome.return_value = "GCA000000000"
        mock_reference_fasta = "assembly.fasta"
        mock_assembly_report = "assembly_report.txt"
        mock_platform.return_value = "myplatform"
        mock_software.return_value = ["software1", "software2"]
        mock_pipeline_descriptions.return_value = "my pipeline description"
        mock_run_accessions.return_value = ["ERR123456", "SRR123456"]
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_analysis("estd123", mock_reference_fasta, mock_assembly_report)
        expected_result = [{'analysisTitle': 'estd123',
                            'analysisAlias': 'MYanalysisALIAS_ALL',
                            'description': 'mock_desc',
                            'experimentType': 'genotyping_by_sequencing',
                            'referenceGenome': 'GCA000000000',
                            'referenceFasta': mock_reference_fasta,
                            'assemblyReport': mock_assembly_report,
                            'platform': 'myplatform',
                            'software': ['software1', 'software2'],
                            'pipelineDescriptions': 'my pipeline description',
                            'runAccessions': ['ERR123456', 'SRR123456']}]
        self.assertEqual(result,expected_result)
        # testing required only
        mock_ids.return_value = ["1","2", "3"]
        mock_alias.return_value = "MYanalysisALIAS_REQUIRED"
        mock_desc.return_value = "mock_desc"
        mock_analysis_type.return_value = "Genotyping"
        mock_method_type.return_value = "Sequencing"
        mock_determined_experiment_type.return_value = "genotyping_by_sequencing"
        mock_reference_genome.return_value = "GCA000000000"
        mock_reference_fasta = "assemble.fasta"
        mock_assembly_report = "assembly_report.txt"
        mock_platform.return_value = ""
        mock_software.return_value = []
        mock_pipeline_descriptions.return_value = ""
        mock_run_accessions.return_value = []
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_analysis("estd123", mock_reference_fasta, mock_assembly_report)

        expected_result = [{'analysisTitle': 'estd123',
                            'analysisAlias': 'MYanalysisALIAS_REQUIRED',
                            'description': 'mock_desc',
                            'experimentType': 'genotyping_by_sequencing',
                            'referenceGenome': 'GCA000000000',
                            'referenceFasta': mock_reference_fasta,
                            'assemblyReport': mock_assembly_report,
        }]
        self.assertEqual(result, expected_result)

    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_sample_analysis_alias_list")
    def test__get_sample_pre_registered(self, mock_alias_list):
        mock_alias_list.return_value = ["alias"]
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_sample_pre_registered("STUDY123", "SAMNA123456", "favourite_sample")
        expected = {'analysisAlias': ['alias'], 'sampleInVCF': 'favourite_sample', 'bioSampleAccession': 'SAMNA123456'}
        self.assertEqual(result, expected)


    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_scientific_name")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_tax_id")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_sample_analysis_alias_list")
    def test__get_sample_new(self, mock_alias_list, mock_tax_id, mock_scientific_name):
        mock_alias_list.return_value = ["alias"]
        mock_tax_id.return_value = '9606'
        mock_scientific_name.return_value = "homo sapiens"

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_sample_new("STUDY123", "favourite_sample")
        expected = {
            'analysisAlias': ['alias'],
            'sampleInVCF': 'favourite_sample',
            'bioSampleObject': {
                'name': 'favourite_sample',
                'characteristics': {
                    'organism': [{'text': 'homo sapiens'}],
                    'taxId': [{'text': '9606'}],
                    'collection date': [{'text': 'not provided'}],
                    'geographic location (country and/or sea)': [{'text': 'not provided'}]
                }
            }
        }
        self.assertEqual(result, expected)


    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_alias")
    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever._fetch_analysis_ids")
    def test__get_files(self, mock_ids, mock_alias):
        mock_ids.return_value = ["1", "2"]
        mock_alias.return_value = "alias_1_2"

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._get_files("STUDY123")
        expected = [{'analysisAlias': 'alias_1_2', 'fileName': ''}] # only required fields
        self.assertEqual(result[0].get("analysisAlias"), expected[0].get("analysisAlias"))
        self.assertEqual(result[0].get("fileName"), expected[0].get("fileName"))


    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever.load_from_db")
    def test__fetch_project_peer_project(self, mock_load):
        # create mock data from load_from_db (PROJECT NOT PRE-REGISTERED)
        mock_data = [(None,) ] # not pre-registered
        mock_load.return_value = mock_data

        metadata_client = EVAMetadataRetriever(self.config)
        accession = metadata_client._fetch_project_peer_project("STUDY123")
        self.assertEqual(accession, None)
        # create mock data from load_from_db (PROJECT PRE-REGISTERED)
        mock_data = [('PRJNA28889',)] # pre-registered

        mock_load.return_value = mock_data
        metadata_client = EVAMetadataRetriever(self.config)
        accession = metadata_client._fetch_project_peer_project("STUDY789")
        self.assertEqual(accession, "PRJNA28889")

    @patch("convert_gvf_to_vcf.metadata_retrievers.eva_metadata.EVAMetadataRetriever.load_from_db")
    def test__determine_sample_pre_registered(self, mock_load):
        # (SAMPLE NOT PRE-REGISTERED)
        # create mock data from load_from_db (SAMPLE NOT PRE-REGISTERED)
        mock_data = [ (None, "NA20344"), (None, "NA20345")  ]# not pre-registered
        mock_load.return_value = mock_data

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._determine_sample_pre_registered("estd123")

        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])
        expected_result = [SampleStatus(is_sample_preregistered=False, sample_accession=None, sample_id="NA20344"), SampleStatus(is_sample_preregistered=False, sample_accession=None, sample_id="NA20345")]
        self.assertEqual(result, expected_result)

        # create mock data from load_from_db (SAMPLE PRE-REGISTERED)
        mock_data = [("SAMN12345678","mypreregisteredsample")]
        mock_load.return_value = mock_data

        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._determine_sample_pre_registered("estd123")

        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])
        expected_result = [SampleStatus(is_sample_preregistered=True, sample_accession="SAMN12345678", sample_id="mypreregisteredsample")]
        self.assertEqual(result, expected_result)

    def test__determine_analysis_experiment_type(self):
        analysis_types= ["Read depth and paired-end mapping"]
        method_types = ["Sequencing"]
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._determine_analysis_experiment_type(analysis_types, method_types)
        expected = "Whole genome sequencing"
        self.assertEqual(result, expected)


    def test__format_project(self):
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client._format_project_fields(None, ["www.link.com", "www.link2.com"], None, "123456")
        expected = ('', ['www.link.com| URL', 'www.link2.com| URL'], '', ['PubMed:123456'])
        self.assertEqual(result, expected)

    def test_validate_fetch_result(self):
        metadata_client = EVAMetadataRetriever(self.config)
        # expect string
        result = metadata_client.fetch_results_from_rows("name", [("value",)])
        expected = ["value"]
        self.assertEqual(result, expected)
        # expect multiple values in list
        result = metadata_client.fetch_results_from_rows("name", [("value",)])
        expected = ["value"]
        self.assertEqual(result, expected)

    def test_validate_date(self):
        metadata_client = EVAMetadataRetriever(self.config)
        result = metadata_client.validate_date("2012-12-12")
        expected = True
        self.assertEqual(result,expected)
        result = metadata_client.validate_date("abcdefg")
        expected = False
        self.assertEqual(result, expected)

    def test_validate_project(self):
        metadata_client = EVAMetadataRetriever(self.config)
        description = "a description"
        hold_date = "2012-12-12"
        parent = "PRJEA123456"
        tax_id = 9606
        title = "title"
        publication = ["PubMed:123345"]
        centre = "A new centre."
        child = "PRJEA223456"
        peer = "PRJEA323456"
        links = [""]
        project_object = {
            "title": title,
            "description": description,
            "taxId": tax_id,
            "centre": centre
        }
        project_object_not_required_all = {
            "publications": publication,
            "parentProject": parent,
            "childProject": child,
            "peerProject": peer,
            "links": links,
            "holdDate": hold_date
        }
        required_result = metadata_client.is_project_valid(project_object)
        self.assertEqual(required_result, True)
        project_object.update(project_object_not_required_all)
        with_not_required_result = metadata_client.is_project_valid(project_object)
        self.assertEqual(with_not_required_result, True)

    def test_validate_analysis(self):
        pass