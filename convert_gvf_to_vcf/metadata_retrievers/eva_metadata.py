import hashlib
import json
import os
import re
from datetime import datetime
from collections import namedtuple

import oracledb
from pypika import Query, Table, Schema
from jsonschema import validate, ValidationError

from convert_gvf_to_vcf.project_paths import ProjectPaths
from convert_gvf_to_vcf.utils import read_in_json_schema
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

logger = log_cfg.get_logger(__name__)

from convert_gvf_to_vcf.metadata_retrievers.base_metadata_retriever import BaseMetadataRetriever
class EVAMetadataRetriever(BaseMetadataRetriever):
    """
    The responsibility of this class is to retrieve the metadata for submission to EVA.
    It will extract the required properties to fill in the EVA JSON submission schema
    (https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json)
    by querying the DGVa database.
    It will generate the metadata submission JSON file.
    """

    def create_json_eva(self, json_file_path, study_accession, assembly, assembly_report):
        project_metadata = self._get_project_new(study_accession)  # (all projects are new projects)
        sample_registration_statuses = self._determine_sample_pre_registered(study_accession)
        sample_metadata_array = []

        for sample_status in sample_registration_statuses:
            if sample_status.is_sample_preregistered:
                sample_metadata = self._get_sample_pre_registered(study_accession, sample_status.sample_accession, sample_status.sample_id)
                sample_metadata_array.append(sample_metadata)
            else:
                # assuming all samples are not pre-registered
                sample_metadata = self._get_sample_new(study_accession, sample_status.sample_id)
                sample_metadata_array.append(sample_metadata)
        json_in_eva_format = {
            "submitterDetails": self._get_submitter_details(study_accession),
            "project": project_metadata,
            "analysis": self._get_analysis(study_accession, assembly, assembly_report),
            "sample": sample_metadata_array,
            "files": self._get_files(study_accession)
        }
        os.makedirs(os.path.dirname(json_file_path), exist_ok=True)
        with open(json_file_path, 'w') as f:
            json.dump(json_in_eva_format, f, indent=4)
        logger.info(f"Write JSON file for {study_accession}- SUCCESS: {json_file_path}")

    # THESE GETTERS GET THE RELEVANT JSON OBJECT FOR THE SECTION
    def _get_submitter_details(self, study_accession):
        """ Formats the submitter details from the DGVa database. Creates a submitter object for each submitter.
        :param: study_accession - expected format ^(estd|nstd)\d+$
        :return: submitter_details_array
        """
        logger.info("Fetching submitter details.")

        all_last_names = self._fetch_submitter_details_all_last_names(study_accession)
        all_first_names = self._fetch_submitter_details_all_first_names(study_accession)
        all_phone_numbers = self._fetch_submitter_details_all_phone_numbers(study_accession)
        all_email_addresses = self._fetch_submitter_details_all_email_addresses(study_accession)
        all_centres = self._fetch_submitter_details_all_centres(study_accession)
        all_addresses = self._fetch_submitter_details_all_addresses(study_accession)
        # load metadata - expected values in form of LIST OF TUPLES [(TUPLE FOR SUBMITTER 1), (TUPLE FOR SUBMITTER 2)]
        # TUPLE = ("FIRST_NAME", "LAST_NAME", "TELEPHONE", "EMAIL", "CENTRE", "ADDRESS")
        # required
        submitter_details_array = [
            {
                "lastName": last_name or "",
                "firstName": first_name or "",
                "email": email_address or "",
                "laboratory": "",  # expected to be empty because this is not found in dgva
                "centre": centre or "",
            }
            for last_name, first_name, email_address, centre in zip(all_last_names, all_first_names, all_email_addresses, all_centres)]

        # not required
        submitter_details_array_not_required_all = [
            {
            "telephone": phone_number,
            "address": address
            } for phone_number, address in zip(all_phone_numbers, all_addresses)]
        submitter_details_array_not_required = []
        for submitter in submitter_details_array_not_required_all:
            not_required = {k: v for k, v in submitter.items() if v}
            submitter_details_array_not_required.append(not_required)
        for required_dict, not_required_dict in zip(submitter_details_array, submitter_details_array_not_required):
            required_dict.update(not_required_dict)

        return submitter_details_array

    def _get_project_pre_registered(self, project_accession):
        # check project accession meets the regex
        project_accession_pattern = r"^PRJ(E|D|N)[A-Z][0-9]+$"
        assert re.fullmatch(project_accession_pattern, project_accession), f"Project accession {project_accession} does not meet the regex: {project_accession_pattern}"
        logger.info(f"Project accession {project_accession} has been found as pre-registered.")
        project_object = {
            "projectAccession" : project_accession
        }
        return project_object

    def _get_project_new(self, study_accession):
        logger.info("Fetching Project details.")
        project_title = self._fetch_project_title(study_accession)
        project_description = self._fetch_project_description(study_accession)
        project_tax_id = self._fetch_tax_id(study_accession)
        project_centre_per_submitter = self._fetch_centre(study_accession)
        project_centre = project_centre_per_submitter[0] if project_centre_per_submitter else "" # choose the first submitter's centre

        project_publications = self._fetch_project_publications(study_accession)
        project_parent_project = self._fetch_project_parent_project(study_accession)
        project_child_project = "" # expect no child projects to exist
        project_peer_project = self._fetch_project_peer_project(study_accession)
        # expect no child projects to exist
        project_links = self._fetch_project_links(study_accession) # expect this to be only one link
        project_hold_date = self._fetch_hold_date(study_accession)

        # formatting
        project_hold_date, project_links, project_parent_project, pubmed_publications = self._format_project_fields(
            project_hold_date, project_links, project_parent_project, project_publications)

        logger.info(f"Project accession has not been found. Creating a new project. ")
        # required: title, description, taxID, centre
        # not required: publications, parentProject, childProject, peerProjects, hold-date, links
        project_object = {
            "title": project_title,
            "description": project_description,
            "taxId": project_tax_id,
            "centre": project_centre
        }
        project_object_not_required_all = {
            "publications": pubmed_publications,
            "parentProject": project_parent_project,
            "childProject": project_child_project,
            "peerProject": project_peer_project,
            "links": project_links,
            "holdDate": project_hold_date
        }
        project_object_not_required = {k:v for k,v in project_object_not_required_all.items() if v}
        project_object.update(project_object_not_required)

        # ensure compliance with EVA JSON schema https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json
        if not self.is_project_valid(project_object):
            raise ValueError("Project does not match the JSON schema")
        return project_object

    def _get_analysis(self, study_accession, assembly, assembly_report):
        # return analysis_array
        # required: analysisTitle, analysisAlias, description, experimentType, reference_genome
        logger.info("Fetching Analysis details.")
        analysis_analysis_id_list = self._fetch_analysis_ids(study_accession)
        # analysis alias is a string in Analysis section
        analysis_analysis_alias = self._fetch_analysis_alias(study_accession, analysis_analysis_id_list)
        analysis_analysis_title = study_accession
        analysis_analysis_description = self._fetch_analysis_description(study_accession)
        analysis_types = self._fetch_analysis_analysis_type(study_accession)
        method_types = self._fetch_analysis_method_type(study_accession)
        analysis_experiment_type = self._determine_analysis_experiment_type(analysis_types, method_types)
        analysis_reference_genome = self._fetch_analysis_reference_genome(study_accession)
        analysis_evidence_type = ""
        analysis_reference_fasta = assembly
        analysis_assembly_report = assembly_report
        analysis_platform = self._fetch_analysis_platform(study_accession)
        analysis_software = self._fetch_analysis_software(study_accession)
        analysis_pipeline_descriptions = self._fetch_analysis_pipeline_descriptions(study_accession)
        analysis_imputation = "" # deliberately empty string to avoid assumption
        analysis_phasing = "" #  deliberately empty string to avoid assumption
        analysis_date = "" # not present in DGVA
        analysis_centre = "" # if different to project centre
        analysis_links = "" # these are placed in project links
        analysis_run_accessions = self._fetch_analysis_run_accessions(study_accession)

        # cleaning and validation
        analysis_pipeline_descriptions, analysis_run_accessions = self.validate_analysis(analysis_pipeline_descriptions,
                                                                                         analysis_run_accessions)

        # required: analysisTitle, analysisAlias, description, experimentType, referenceGenome
        # not required: evidenceType, referenceFasta, assemblyReport, platform, software, pipelineDescriptions
        # imputation, phasing, date, links, runAccessions
        analysis_array = []
        analysis_object = {
            "analysisTitle": analysis_analysis_title,
            "analysisAlias": analysis_analysis_alias,
            "description": analysis_analysis_description,
            "experimentType": analysis_experiment_type,
            "referenceGenome": analysis_reference_genome,
        }
        analysis_object_not_required_all = {
            "evidenceType": analysis_evidence_type,
            "referenceFasta": analysis_reference_fasta,
            "assemblyReport": analysis_assembly_report,
            "platform": analysis_platform,
            "software": analysis_software,
            "pipelineDescriptions": analysis_pipeline_descriptions,
            "imputation": analysis_imputation,
            "phasing": analysis_phasing,
            "date": analysis_date,
            "centre": analysis_centre,
            "links":analysis_links,
            "runAccessions": analysis_run_accessions
        }
        analysis_object_not_required = {k: v for k, v in analysis_object_not_required_all.items() if v}
        analysis_object.update(analysis_object_not_required)
        analysis_array.append(analysis_object)
        return analysis_array

    def _get_sample_pre_registered(self, study_accession, biosample_accession, sample_id):
        # requires analysisAlias, sampleinVCF, biosample_accession
        # analysis_alias is an array in samples
        sample_analysis_alias_list = self._fetch_sample_analysis_alias_list(study_accession, sample_id)
        # assumption the name: sampleinVCF = sample_id
        sample_sampleinvcf = sample_id
        sample_object = {
            "analysisAlias": sample_analysis_alias_list,
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleAccession": biosample_accession
        }
        return sample_object

    def _get_sample_new(self, study_accession, sample_id):
        # return sample object
        # requires analysisAlias, sampleinVCF, bioSampleObject
        # bioSampleObject requires = name, taxID, scientific_name, release (hold-date) which can be found in DGVA
        # bioSampleObject requires = collection date, geo loc which can be set to unknown/not collected
        sample_analysis_alias_list = self._fetch_sample_analysis_alias_list(study_accession, sample_id)
        if not sample_analysis_alias_list:
            sample_analysis_alias_list.append("")
        # assume sample in VCF = sample_id
        sample_sampleinvcf = sample_id
        sample_tax_id = str(self._fetch_tax_id(study_accession))
        scientific_name = self._fetch_scientific_name(study_accession)
        # hold_date = self._fetch_hold_date(study_accession)
        collection_date = "not provided"
        geographic_location_country_and_or_sea = "not provided"
        biosample_object = {
            "name": sample_id,
            "characteristics": {
                "organism": [{"text": scientific_name}],
                "taxId":[{"text": sample_tax_id}],
                "collection date": [{"text": collection_date}],
                "geographic location (country and/or sea)": [{"text": geographic_location_country_and_or_sea}]
            }
        }
        sample_object = {
            "analysisAlias": sample_analysis_alias_list,
            "sampleInVCF": sample_sampleinvcf,
            "bioSampleObject": biosample_object
        }
        return sample_object

    def _get_files(self, study_accession):
        files_analysis_id_list = self._fetch_analysis_ids(study_accession)
        # analysis alias is a string in files
        files_analysis_alias = self._fetch_analysis_alias(study_accession, files_analysis_id_list)

        files_file_name = ""
        files_file_size = ""
        files_file_md5 = ""

        files_array = []
        # required: analysisAlias, filename
        # not required: file size, md5
        files_object = {
            "analysisAlias": files_analysis_alias,
            "fileName": files_file_name
        }
        files_object_not_required_all = {
            "fileSize": files_file_size,
            "md5": files_file_md5
        }
        files_object_not_required = {k: v for k, v in files_object_not_required_all.items() if v}
        files_object.update(files_object_not_required)
        files_array.append(files_object)
        return files_array

    # THESE DETERMINE IF NEW OR PRE_REG OR EVIDENCE_TYPE OR EXPERIMENT_TYPE


    def _determine_sample_pre_registered(self, study_accession):
        # create the table objects
        dsamp = Table("DGVA_SAMPLE", schema=self.db).as_("dsamp")
        # create the query
        sample_accession_query = (Query
                                   .from_(dsamp)
                                   .select(dsamp.BIOSAMPLE_ACCESSION, dsamp.SUBMITTER_SAMPLE_ID)
                                   .where(dsamp.STUDY_ACCESSION == study_accession)
                                   )
        # sample_accession_list: [(BiosampleAccession, SubmitterSampleID)]
        sample_accession_list = self.load_from_db(sample_accession_query.get_sql(quote_char=None))
        sample_accession_and_status_list = []
        SampleStatus = namedtuple('SampleStatus', ['is_sample_preregistered', 'sample_accession', 'sample_id'])

        is_sample_preregistered = False
        for value in sample_accession_list:
            if isinstance(value, tuple) and len(value) > 0:
                current_sample_accession = value[0]
                current_sample_id = value[1]
                if current_sample_accession is not None:
                    # DGVa does not have any samples with a BioSample accession. We do not expect this to be populated.
                    is_sample_preregistered = True
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, current_sample_accession, current_sample_id))
                    # return is_sample_preregistered, sample_accession, sample_id
                    logger.info(f"Determining if sample is pre-registered - SUCCESS - Sample found: {SampleStatus(is_sample_preregistered, current_sample_accession, current_sample_id)}.")
                else:
                    sample_accession_and_status_list.append(SampleStatus(is_sample_preregistered, None, current_sample_id))
                    # return is_sample_preregistered, None, sample_id
                    logger.info(f"Determining if sample is pre-registered - FAILURE - Sample not found: {current_sample_id}.")
        return sample_accession_and_status_list

    def _determine_analysis_experiment_type(self, analysis_types, method_types):
        analysis_and_method_types = list(zip(analysis_types, method_types))
        # mapping method type (value) to experiment type (key)
        experiment_type_to_method_type = {
            "whole_genome_sequencing": ["Sequencing", "Merging"],
            "whole_transcriptome_sequencing": ["Sequencing", "Merging"],
            "exome_sequencing": ["Sequencing", "Merging"],
            "genotyping_by_array": ["SNP array", "Oligo aCGH", "MLPA", "FISH", "Karyotyping"],
            "genotyping_by_sequencing": ["Sequencing", "PCR", "qPCR", "Digital array"],
            "target_sequencing": ["Sequencing", "PCR", "qPCR", "MLPA"],
            "transcriptomics": ["Sequencing", "qPCR", "Digital array"],
            "curation": ["Merging"]
        }
        # mapping analysis type (key) to experiment type (value)
        analysis_to_experiment = {
            "Read depth and paired-end mapping": "Whole genome sequencing",
            "Manual observation": "Curation",
            "de novo sequence assembly": "Whole genome sequencing",
            "Sequence alignment": "Whole genome sequencing",
            "Genotyping": "Genotyping by sequencing",
            "Split read mapping": "Whole genome sequencing",
            "Paired-end mapping": "Whole genome sequencing",
            "Split read and paired-end mapping": "Whole genome sequencing",
            "Merging": "Curation",
            "Other": "Curation",
            "SNP genotyping analysis": "Genotyping by array",
            "Probe signal intensity": "Genotyping by array"
        }
        # key = experiment type; value = experiment type formatted for schema
        experiment_type_to_schema_experiment_type = {
            "whole_genome_sequencing": "Whole genome sequencing",
            "whole_transcriptome_sequencing": "Whole transcriptome sequencing",
            "exome_sequencing": "Exome sequencing",
            "genotyping_by_array": "Genotyping by array",
            "genotyping_by_sequencing": "Genotyping by sequencing",
            "target_sequencing": "Target sequencing",
            "transcriptomics": "Transcriptomics",
            "curation": "Curation"
        }
        experiment_type_set = set()
        for i in analysis_and_method_types:
            analysis_type = i[0]
            method_type = i[1]
            raw_exp_type = analysis_to_experiment.get(analysis_type)
            exp_type = (raw_exp_type or "").lower().replace(" ", "_")
            accepted_method_types = experiment_type_to_method_type.get(exp_type, [])
            if method_type in accepted_method_types:
                experiment_type_set.add(exp_type)
        if len(experiment_type_set) == 1:
            experiment_type = list(experiment_type_set)[0]
            return experiment_type_to_schema_experiment_type.get(experiment_type)
        return None

    # Formatting
    def create_set(self, list_of_elements):
        set_of_elements = set()
        for element in (list_of_elements or []):
            if element is not None:
                set_of_elements.add(element)
        return set_of_elements

    def concatenate_elements_to_string(self, set_or_list):
        string_to_generate = ""
        for element in set_or_list:
            string_to_generate += " & " + str(element)
        string_to_generate = string_to_generate.lstrip(" & ")
        return string_to_generate

    def _format_project_fields(self, project_hold_date, project_links, project_parent_project, project_publications):
        if project_hold_date is None:
            project_hold_date = ""
        label = "| URL"
        if project_links != "" and project_links is not None:
            for index, link in enumerate(project_links):
                project_links[index] =  f"{link}{label}" if link else ""
                # project_links = project_links + "| URL"
        else:
            project_links = ""
        pubmed_publications = []
        if project_parent_project is None:
            project_parent_project = ""
        if type(project_publications) is not list:
            project_publications = [project_publications]
        for pub in project_publications:
            if pub:
                pubmed_string = "PubMed:" + str(pub)
                pubmed_publications.append(pubmed_string)
        return project_hold_date, project_links, project_parent_project, pubmed_publications

    def validate_date(self, date):
        try:
            datetime.strptime(date, '%Y-%m-%d')
            return True
        except ValueError:
            return False

    def is_project_valid(self, project_object_to_validate):
        """ Validates project with the EVA JSON schema https://github.com/EBIvariation/eva-sub-cli/blob/main/eva_sub_cli/etc/eva_schema.json
        :params: project_description: string of max 5000 chars
        :params: project_hold_date: YYYY-MM-DD or ""
        :params: project_parent_project: project accession matching regex "^PRJ(E|D|N)[A-Z][0-9]+$"
        :params: project_tax_id: taxonomy ID as an integer
        :params: project_title: string of max 500 chars
        :params: pubmed_publications: list of pubmed ids ['Pubmed:1239234'] (can be more than one, in some studies)
        """
        schema = read_in_json_schema(ProjectPaths().eva_schema_path)
        project_schema = schema["properties"]["project"]
        project_schema["definitions"] = schema["definitions"]
        try:
            validate(instance=project_object_to_validate, schema=project_schema)
            return True
        except (FileNotFoundError, json.JSONDecodeError, ValidationError) as e:
            logger.error(f"Validating Project Error: {e}")
            return False

    def validate_analysis(self, analysis_pipeline_descriptions, analysis_run_accessions):
        if analysis_run_accessions is not None:
            analysis_run_accessions = [r if r else "" for r in analysis_run_accessions]
            run_accession_pattern = "^(E|D|S)RR[0-9]{6,}$"
            for run_accession in analysis_run_accessions:
                if run_accession is not None and run_accession != "":
                    assert re.fullmatch(run_accession_pattern,
                                        run_accession), f"String {run_accession} does not match pattern: {run_accession_pattern}"
            if all(not run_accession for run_accession in analysis_run_accessions):
                analysis_run_accessions = ""
        if all(not pipeline_descriptions for pipeline_descriptions in analysis_pipeline_descriptions):
            analysis_pipeline_descriptions = ""
        return analysis_pipeline_descriptions, analysis_run_accessions

    #_fetch_<section>_all_fields = return a list
    #_fetch_<section>_field_name = return a single value
    # SUBMITTER DETAILS SECTION
    def _fetch_submitter_details_all_last_names(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_last_names_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.LAST_NAME
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_last_names_list = self.load_from_db(all_last_names_query.get_sql(quote_char=None))
        all_last_names = self.fetch_results_from_rows("lastName", all_last_names_list)
        return all_last_names

    def _fetch_submitter_details_all_first_names(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_first_names_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.FIRST_NAME
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_first_names_list = self.load_from_db(all_first_names_query.get_sql(quote_char=None))
        all_first_names = self.fetch_results_from_rows("firstName", all_first_names_list)
        return all_first_names

    def _fetch_submitter_details_all_phone_numbers(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_phone_numbers_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.CONTACT_PHONE
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_phone_numbers_list = self.load_from_db(all_phone_numbers_query.get_sql(quote_char=None))
        all_phone_numbers = self.fetch_results_from_rows("telephone", all_phone_numbers_list)
        return all_phone_numbers

    def _fetch_submitter_details_all_email_addresses(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_email_addresses_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.CONTACT_EMAIL
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_email_addresses_list = self.load_from_db(all_email_addresses_query.get_sql(quote_char=None))
        all_email_addresses = self.fetch_results_from_rows("email", all_email_addresses_list)
        return all_email_addresses

    def _fetch_submitter_details_all_centres(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_centres_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.AFFILIATION_NAME
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_centres_list = self.load_from_db(all_centres_query.get_sql(quote_char=None))
        all_centres = self.fetch_results_from_rows("centre", all_centres_list)
        return all_centres

    def _fetch_submitter_details_all_addresses(self, study_accession):
        # create the table objects
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # programmatically create the queries
        all_addresses_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(
                sc.AFFILIATION_ADDRESS
            ).where(ds.STUDY_ACCESSION == study_accession)
        )
        all_addresses_list = self.load_from_db(all_addresses_query.get_sql(quote_char=None))
        all_addresses = self.fetch_results_from_rows("address", all_addresses_list)
        return all_addresses
    # PROJECT SECTION
    def _fetch_project_title(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        project_title_query = (
            Query.from_(ds)
            .select(ds.DISPLAY_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)

        )
        project_title_list = self.load_from_db(project_title_query.get_sql(quote_char=None))
        [project_title, *_] = self.fetch_results_from_rows("title", project_title_list) or [""]
        return project_title

    def _fetch_project_description(self, study_accession):
        # create the table objects
        stc = Table("STUDY_TYPE_CV", schema=self.db).as_("stc")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        project_description_query = (
            Query.from_(ds)
            .join(stc)
            .on(ds.STUDY_TYPE == stc.STUDY_TYPE)
            .select(ds.STUDY_DESCRIPTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        project_description_list = self.load_from_db(project_description_query.get_sql(quote_char=None))
        [project_description, *_] = self.fetch_results_from_rows("description", project_description_list) or [""]
        return project_description

    def _fetch_tax_id(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        so = Table("STUDY_ORGANISM", schema=self.db).as_("so")
        tax_id_query = (
            Query.from_(ds)
            .join(so)
            .on(so.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(so.TAXONOMY_ID)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        tax_id_list = self.load_from_db(tax_id_query.get_sql(quote_char=None))
        [tax_id, *_] = self.fetch_results_from_rows("taxId", tax_id_list) or [""]
        return tax_id

    def _fetch_centre(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sc = Table("STUDY_CONTACT", schema=self.db).as_("sc")
        project_centre_query = (
            Query.from_(sc)
            .join(ds)
            .on(sc.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(sc.AFFILIATION_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        project_centre_list = self.load_from_db(project_centre_query.get_sql(quote_char=None))
        [project_centre, *_] = self.fetch_results_from_rows("centre", project_centre_list) or [""]
        return project_centre

    def _fetch_project_publications(self, study_accession):
        # some study accessions do have multiple publications e.g. estd192
        # create the table objects
        spp = Table("STUDY_PUBMED_PUBLICATION", schema=self.db).as_("spp")
        project_publications_query = (
            Query.from_(spp)
            .select(spp.PUBMED_ID)
            .where(spp.STUDY_ACCESSION == study_accession)
        )
        project_publications_list = self.load_from_db(project_publications_query.get_sql(quote_char=None))
        [project_publications, *_] = self.fetch_results_from_rows("publications", project_publications_list) or [""]
        return project_publications

    def _fetch_project_parent_project(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # create the query
        parent_project_query = (Query
                                   .from_(ds)
                                   .select(ds.BIOPROJECT_ACCESSION)
                                   .where(ds.STUDY_ACCESSION == study_accession)
                                   )
        parent_project_list = self.load_from_db(parent_project_query.get_sql(quote_char=None))
        [parent_project, *_] = self.fetch_results_from_rows("projectAccession", parent_project_list) or [""]
        return parent_project

    def _fetch_project_peer_project(self, study_accession):
        """ Determines if the project is new or pre-registered.
        :param: study_accession - expected format ^(estd|nstd)\d+$
        :return project_accession : string format ^(PRJ)[A-Z]{2}\d+$
        """
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        # create the query
        peer_project_accession_query = (Query
                                   .from_(ds)
                                   .select(ds.BIOPROJECT_ACCESSION)
                                   .where(ds.STUDY_ACCESSION == study_accession)
                                   )
        peer_project_accession_list = self.load_from_db(peer_project_accession_query.get_sql(quote_char=None))
        [peer_project_accession, *_] = self.fetch_results_from_rows("peerProject", peer_project_accession_list) or [""]
        return peer_project_accession

    def _fetch_project_links(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        sub = Table("DGVA_SUBMISSION", schema=self.db).as_("sub")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        elu = Table("EXPERIMENT_LINK_URL", schema=self.db).as_("elu")
        ec = Table("EXPERIMENT_CURATION", schema=self.db).as_("ec")
        # TABLE_NAME:COLUMN_NAME = dgva_study:study_url
        project_links_query_one = (
            Query.from_(ds)
            .join(sub)
            .on(ds.STUDY_ACCESSION == sub.STUDY_ACCESSION)
            .select(ds.STUDY_URL)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        # TABLE_NAME:COLUMN_NAME = experiment_link_url:url
        project_links_query_two = (
            Query.from_(elu)
            .join(de)
            .on(de.EXPERIMENT_ID == elu.EXPERIMENT_ID)
            .select(elu.URL)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        # TABLE_NAME:COLUMN_NAME = experiment_curation:curated_set_link
        project_links_query_three = (
            Query.from_(ec)
            .join(de)
            .on(de.EXPERIMENT_ID == ec.EXPERIMENT_ID)
            .select(ec.CURATED_SET_LINK)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        project_links_queries = [project_links_query_one, project_links_query_two, project_links_query_three]
        project_links = set()
        for query in project_links_queries:
            query_links_list = self.load_from_db(query.get_sql(quote_char=None))
            fetched_results = self.fetch_results_from_rows("links", query_links_list)
            if fetched_results:
                project_links.update(fetched_results)

        # adding analysis links here because of project_links flexibility
        analysis_links = self._fetch_analysis_links(study_accession)
        if analysis_links:
            project_links.update(analysis_links)
        project_links_list = list(project_links)
        return project_links_list

    # ANALYSIS SECTION
    def _fetch_analysis_ids(self, study_accession):
        # create the table objects
        da = Table("DGVA_ANALYSIS", schema=self.db).as_("da")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ea = Table("EXPERIMENT_ANALYSIS", schema=self.db).as_("ea")
        analysis_id_query = (
            Query.from_(de)
            .select(da.ANALYSIS_ID)
            .join(ea).on(ea.EXPERIMENT_ID == de.EXPERIMENT_ID )
            .join(da).on(da.ANALYSIS_ID == ea.ANALYSIS_ID)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_id_list = self.load_from_db(analysis_id_query.get_sql(quote_char=None))
        all_analysis_ids = self.fetch_results_from_rows("analysisID", analysis_id_list)
        all_analysis_ids.sort()
        return all_analysis_ids

    def _fetch_analysis_alias(self, study_accession, analysis_id_list):
        if len(analysis_id_list) == 1:
            analysis_prefix = "_analysis_id_"
        else:
            analysis_prefix = "_analysis_ids_"
        analysis_alias = study_accession + analysis_prefix
        for analysis_id in analysis_id_list:
            analysis_alias = analysis_alias + str(analysis_id) + "_"
        analysis_alias = analysis_alias.rstrip("_")
        return analysis_alias

    def _fetch_analysis_description(self, study_accession):
        # create the table objects
        da = Table("DGVA_ANALYSIS", schema=self.db).as_("da")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ea = Table("EXPERIMENT_ANALYSIS", schema=self.db).as_("ea")

        analysis_description_query = (
            Query.from_(ds)
            .join(de).on(de.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .join(ea).on(de.EXPERIMENT_ID == ea.EXPERIMENT_ID)
            .join(da).on(ea.ANALYSIS_ID == da.ANALYSIS_ID)
            .select(da.ANALYSIS_DESCRIPTION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        analysis_description_list = self.load_from_db(analysis_description_query.get_sql(quote_char=None))
        analysis_descriptions = self.fetch_results_from_rows("description", analysis_description_list)
        analysis_description = self.concatenate_elements_to_string(analysis_descriptions)
        return analysis_description

    def _fetch_analysis_analysis_type(self, study_accession):
        # create the table objects
        vse = Table("VW_SUMM_EXPERIMENT", schema=self.db).as_("vse")
        analysis_type_query = (
            Query.from_(vse)
            .select(vse.ANALYSIS_TYPE)
            .where(vse.STUDY_ACCESSION == study_accession)
        )
        analysis_type_list = self.load_from_db(analysis_type_query.get_sql(quote_char=None))
        analysis_types = self.fetch_results_from_rows("analysis_type", analysis_type_list)
        return analysis_types

    def _fetch_analysis_method_type(self, study_accession):
        # create the table objects
        vse = Table("VW_SUMM_EXPERIMENT", schema=self.db).as_("vse")
        analysis_type_query = (
            Query.from_(vse)
            .select(vse.METHOD_TYPE)
            .where(vse.STUDY_ACCESSION == study_accession)
        )
        method_type_list = self.load_from_db(analysis_type_query.get_sql(quote_char=None))
        method_types = self.fetch_results_from_rows("method_types", method_type_list)
        return method_types

    def _fetch_analysis_experiment_type(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        experiment_type_query = (
            Query.from_(de)
            .join(ds).on(de.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(de.EXPERIMENT_TYPE)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        experiment_type_list = self.load_from_db(experiment_type_query.get_sql(quote_char=None))
        experiment_types = self.fetch_results_from_rows("experimentType", experiment_type_list)
        if experiment_types and all(experiment_type == experiment_types[0] for experiment_type in experiment_types):
            experiment_type = experiment_types[0]
            return experiment_type
        return None

    def _fetch_analysis_reference_genome(self, study_accession):
        # create the table objects
        sap = Table("STUDY_ASSEMBLY_PLACEMENT", schema=self.db).as_("sap")
        reference_genome_query = (
            Query.from_(sap)
            .select(sap.ASSEMBLY_ACCESSION)
            .distinct()
            .where(
                (sap.STUDY_ACCESSION == study_accession) & (sap.VARIANT_PLACEMENT_METHOD.isnull())
            )
        )
        reference_genome_list = self.load_from_db(reference_genome_query.get_sql(quote_char=None))
        variant_placement_methods = ["Submitted genomic", "Submitted cytogenic", "Submitted artifact", "Remapped"]
        # if empty, try another variant placement method.
        for variant_placement_method in variant_placement_methods:
            if reference_genome_list == []:
                reference_genome_query = (
                    Query.from_(sap)
                    .select(sap.ASSEMBLY_ACCESSION)
                    .distinct()
                    .where(
                        (sap.STUDY_ACCESSION == study_accession) & (sap.VARIANT_PLACEMENT_METHOD == variant_placement_method)
                    )
                )
                reference_genome_list = self.load_from_db(reference_genome_query.get_sql(quote_char=None))
            if reference_genome_list:
                break
        [reference_genome, *_] = self.fetch_results_from_rows("referenceGenome", reference_genome_list) or [""]
        return reference_genome

    def _fetch_analysis_platform(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ep = Table("EXPERIMENT_PLATFORM", schema=self.db).as_("ep")
        dp = Table("DGVA_PLATFORM", schema=self.db).as_("dp")
        analysis_platform_query = (
            Query.from_(de)
            .join(ep).on(de.EXPERIMENT_ID == ep.EXPERIMENT_ID)
            .join(dp).on(dp.PLATFORM_ID == ep.PLATFORM_ID)
            .select(dp.PLATFORM_NAME)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_platform_list = self.load_from_db(analysis_platform_query.get_sql(quote_char=None))
        analysis_platforms = self.fetch_results_from_rows("platform", analysis_platform_list)
        analysis_platform_set = self.create_set(analysis_platforms)
        analysis_platform = self.concatenate_elements_to_string(analysis_platform_set)
        return analysis_platform

    def _fetch_analysis_software(self, study_accession):
        # create the table objects
        dd = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        analysis_software_query = (
            Query.from_(dd)
            .join(ed).on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de).on(ed.EXPERIMENT_ID == de.EXPERIMENT_ID)
            .select(dd.DETECTION_METHOD)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_software_list = self.load_from_db(analysis_software_query.get_sql(quote_char=None))
        analysis_software = self.fetch_results_from_rows("software", analysis_software_list)
        return analysis_software

    def _fetch_analysis_pipeline_descriptions(self, study_accession):
        # create the table objects
        dd = Table("DGVA_DETECTION", schema=self.db).as_("dd")
        ed = Table("EXPERIMENT_DETECTION", schema=self.db).as_("ed")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        analysis_pipeline_descriptions_query = (
            Query.from_(dd)
            .join(ed).on(dd.DETECTION_ID == ed.DETECTION_ID)
            .join(de).on(ed.EXPERIMENT_ID == de.EXPERIMENT_ID)
            .select(dd.DETECTION_DESCRIPTION)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_pipeline_descriptions_list = self.load_from_db(analysis_pipeline_descriptions_query.get_sql(quote_char=None))
        [analysis_pipeline_descriptions, *_] = self.fetch_results_from_rows("pipelineDescriptions",
                                                                            analysis_pipeline_descriptions_list) or [""]
        analysis_pipeline_description_set = self.create_set(analysis_pipeline_descriptions)
        analysis_pipeline_description = self.concatenate_elements_to_string(analysis_pipeline_description_set)
        return analysis_pipeline_description

    def _fetch_analysis_links(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        elu = Table("EXPERIMENT_LINK_URL", schema=self.db).as_("elu")
        analysis_links_query = (
            Query.from_(de)
            .join(elu).on(de.EXPERIMENT_ID == elu.EXPERIMENT_ID)
            .select(elu.URL)
            .where(de.STUDY_ACCESSION == study_accession)
        )
        analysis_links_list = self.load_from_db(analysis_links_query.get_sql(quote_char=None))
        # NOTE: ANALYSIS_LINKS SHOULD  match regex (DB:ID:LABEL) but we are placing in project links for flexibility.
        analysis_links = self.fetch_results_from_rows("links", analysis_links_list)
        return analysis_links

    def _fetch_analysis_run_accessions(self, study_accession):
        # create the table objects
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        analysis_run_accessions_query = (
            Query.from_(de)
            .join(ds).on(de.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .select(de.SRA_ACCESSION)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        analysis_run_accessions_list = self.load_from_db(analysis_run_accessions_query.get_sql(quote_char=None))
        [analysis_run_accessions, *_] = self.fetch_results_from_rows("runAccessions", analysis_run_accessions_list) or [""]
        return analysis_run_accessions
    # SAMPLE SECTION
    def _fetch_analysis_ids_for_sample(self, study_accession, sample_id):
        # create the table objects
        da = Table("DGVA_ANALYSIS", schema=self.db).as_("da")
        de = Table("DGVA_EXPERIMENT", schema=self.db).as_("de")
        ds = Table("DGVA_SAMPLE", schema=self.db).as_("ds")
        ea = Table("EXPERIMENT_ANALYSIS", schema=self.db).as_("ea")
        analysis_id_for_sample_query = (
            Query.from_(de)
            .select(ea.ANALYSIS_ID)
            .join(ea).on(ea.EXPERIMENT_ID == de.EXPERIMENT_ID )
            .join(da).on(da.ANALYSIS_ID == ea.ANALYSIS_ID)
            .join(ds).on(ds.STUDY_ACCESSION == de.STUDY_ACCESSION)
            .where((de.STUDY_ACCESSION == study_accession) & (ds.SUBMITTER_SAMPLE_ID == sample_id))
        )
        analysis_id_list = self.load_from_db(analysis_id_for_sample_query.get_sql(quote_char=None))
        all_analysis_ids_for_sample = self.fetch_results_from_rows("analysisID", analysis_id_list)
        all_analysis_ids_for_sample.sort()
        return all_analysis_ids_for_sample

    def _fetch_sample_analysis_alias_list(self, study_accession, sample_id):
        analysis_ids_for_sample = self._fetch_analysis_ids_for_sample(study_accession, sample_id)
        # analysis alias is a ARRAY in SAMPLE section
        analysis_alias_list = []
        analysis_analysis_alias_for_sample = self._fetch_analysis_alias(study_accession, analysis_ids_for_sample)
        analysis_alias_list.append(analysis_analysis_alias_for_sample)
        return analysis_alias_list

    def _fetch_sample_id_list(self, study_accession):
        # create the table objects
        dsamp = Table("DGVA_SAMPLE", schema=self.db).as_("dsamp")
        sample_id_query = (
            Query.from_(dsamp)
            .select(dsamp.SUBMITTER_SAMPLE_ID)
            .where(dsamp.STUDY_ACCESSION == study_accession)
        )
        sample_id_list = self.load_from_db(sample_id_query.get_sql(quote_char=None))
        sample_id_keylist = [v[0] for v in sample_id_list]
        return sample_id_keylist

    def _fetch_hold_date(self, study_accession):
        # create the table objects
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        hold_date_query = (
            Query.from_(ds)
            .select(ds.HOLD_DATE)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        hold_date_list = self.load_from_db(hold_date_query.get_sql(quote_char=None))
        [hold_date, *_] = self.fetch_results_from_rows("holdDate", hold_date_list) or [""]
        return hold_date

    def _fetch_scientific_name(self, study_accession):
        # create the table objects
        so = Table("STUDY_ORGANISM", schema=self.db).as_("so")
        ds = Table("DGVA_STUDY", schema=self.db).as_("ds")
        oc = Table("ORGANISM_CV", schema=self.db).as_("oc")
        scientific_name_query = (
            Query.from_(ds)
            .join(so).on(so.STUDY_ACCESSION == ds.STUDY_ACCESSION)
            .join(oc).on(so.TAXONOMY_ID == oc.TAXONOMY_ID)
            .select(oc.SPECIES_LATIN_NAME)
            .where(ds.STUDY_ACCESSION == study_accession)
        )
        scientific_name_list = self.load_from_db(scientific_name_query.get_sql(quote_char=None))
        [scientific_name, *_] = self.fetch_results_from_rows("scientific_name", scientific_name_list) or [""]
        return scientific_name

    # FILES SECTION
    # THESE GETTERS GET THE RELEVANT VALUE
    def _get_file_name(self, vcf_output):
        file_name = os.path.basename(vcf_output)
        study_name = file_name.split(".")[0]
        file_name = os.path.join(study_name, file_name)
        return file_name

    def _get_file_size(self, vcf_output):
        file_size = os.path.getsize(vcf_output)
        return file_size

    def _get_file_md5(self, vcf_output):
        hash_md5 = hashlib.md5()
        with open(vcf_output, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
