import copy
import os
import json
from collections import defaultdict


from convert_gvf_to_vcf.gather_metadata import gather_metadata_workflow, eva_update_metadata_with_vcf
from convert_gvf_to_vcf.convert_gvf_to_vcf_logic import convert, ProjectPaths
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
logger = log_cfg.get_logger(__name__)

class GvfMetadataCoordinator:
    # the responsibility of this class is to co-ordinate multiple GVFs with their Metadata. Feeds into gather metadata script
    def __init__(self, study_and_gvf_files, base_output_dir, config_path):
        self.scan_results = study_and_gvf_files #{study_accession:[gvf_files]}
        self.base_output_dir = base_output_dir # top level output
        self.config_path = config_path
        self.project_paths = ProjectPaths()

    def process_studies(self):
        """ Processes the dictionary result of scanning the directories {study_accession:[gvf_files]}.
        If multiple files it will reconfigure the JSON file separated out by assembly.
        """
        empty_studies_log = []
        for study_accession, gvf_files in self.scan_results.items():
            if gvf_files:
                study_name, _, _ = self.parse_gvf_filename(gvf_files[0])
            else:
                study_name = study_accession
            study_master_json = os.path.join(self.base_output_dir, "submission", study_name,
                                       f"eva_submission_{study_accession}.json")
            # empty list of GVF files
            if not gvf_files:
                self._process_no_gvf_files(empty_studies_log, study_accession)
                continue
            # one file in list of GVFs
            if len(gvf_files) >= 1:
                self._process_gvf_files(gvf_files, study_accession, study_master_json)
            else:
                logger.info("Not in a recognised format")

        logger.info(f"Number of empty studies encountered: {len(empty_studies_log)}")

    @staticmethod
    def _process_no_gvf_files(empty_studies_log, study_accession):
        """Process if no gvf files are present - logging.
        :params: empty_studies_log = list of empty studies
        :params: study_accession e.g. estd1
        """
        logger.info(f"No GVF files found for {study_accession}")
        empty_studies_log.append(study_accession)

    def _process_gvf_files(self, gvf_files, study_accession, master_json):
        """Process if gvf files are present.
        :params: gvf_files = list of gvf files
        :params: study_accession e.g. estd1
        :params: master_json = path to master metadata JSON
        """
        logger.info(f"GVF file(s) found for {study_accession}. Separating by assembly.")
        # separate files by assembly {"GRCh37": [gvf_file_paths], "GRCh38": [gvf_file_paths]})
        assembly_to_gvf_file_paths, gvf_name_groups = self._group_files_by_assembly(gvf_files)
        # master metadata starts off empty and eva_metadata is added to this
        master_metadata = {
            "submitterDetails": None,
            "project": None,
            "analysis": [],
            "sample": [],
            "files": []
        }
        # loop through each assembly
        for assembly_name, gvf_files_for_assembly in assembly_to_gvf_file_paths.items():
            # go through each assembly and process them individually to obtain the metadata
            assembly_path, eva_retriever, json_eva, remapped_files, submitted_files = self._process_single_assembly(
                assembly_name, gvf_files_for_assembly, gvf_name_groups, study_accession
            )
            eva_metadata = self._load_json_file(json_eva)

            if eva_metadata:
                # suffix metadata to prevent same names
                self._suffix_assembly_metadata(eva_metadata, assembly_name)
                # populate master metadata with existing metadata (first metadata file)
                if master_metadata["submitterDetails"] is None and eva_metadata:
                    master_metadata["submitterDetails"] = eva_metadata.get("submitterDetails")
                    master_metadata["project"] = eva_metadata.get("project")
                # overwrite with eva_metadata
                with open(json_eva, 'w', encoding='utf-8') as f_out:
                    json.dump(eva_metadata, f_out)

            for individual_gvf in gvf_files_for_assembly:
                vcf_path = self.convert_individual_gvf(assembly_path, eva_retriever, individual_gvf, json_eva)

            # reconfigure with the master_metadata
            self._reconfigure_assembly_metadata(
                json_eva=json_eva,
                gvf_files=remapped_files if remapped_files else gvf_files_for_assembly,
                master_metadata=master_metadata,
                vcf_path=vcf_path
            )


        if master_metadata["submitterDetails"] is not None:
            self._finalize_master_metadata_cleanup(master_metadata)
            os.makedirs(os.path.dirname(master_json), exist_ok=True)
            with open(master_json, 'w', encoding='utf-8') as master_out:
                json.dump(master_metadata, master_out)
            logger.info(f"Saved aggregated master metadata to {master_json}")

    def _load_json_file(self, file_path):
        """Load a metadata JSON file"""
        if not file_path or not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
            return {}
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError:
            return {}

    def _suffix_assembly_metadata(self, metadata, assembly_name):
        """Add assembly name to the end (analysisAlias, analysisTitle)
        :param metadata = metadata to be changed
        :param assembly_name = name to add
        """
        suffix = f"_{assembly_name}"
        for analysis in metadata.get("analysis", []):
            analysis["analysisAlias"] = f"{analysis.get('analysisAlias', '')}{suffix}"
            analysis["analysisTitle"] = f"{analysis.get('analysisTitle', '')}{suffix}"
        for files in metadata.get("files", []):
            files["analysisAlias"] = f"{files.get('analysisAlias', '')}{suffix}"
        for sample in metadata.get("sample", []):
            sample["analysisAlias"] = [f"{alias}{suffix}" for alias in sample.get("analysisAlias", [])]


    def _finalize_master_metadata_cleanup(self, master_metadata):
        """check for duplicates and merges aliases"""
        if "submitterDetails" in master_metadata and master_metadata["submitterDetails"]:
            seen_emails = set()
            unique_submitters = []
            for submitter in master_metadata["submitterDetails"]:
                email_lower = submitter.get("email", "").lower()
                if email_lower not in seen_emails:
                    seen_emails.add(email_lower)
                    unique_submitters.append(submitter)
            master_metadata["submitterDetails"] = unique_submitters
        merged_samples = {}
        for sample in master_metadata["sample"]:
            sample_name = sample.get("sampleInVCF")
            if sample_name not in merged_samples:
                merged_samples[sample_name] = sample
            else:
                combined_aliases = set(
                    merged_samples[sample_name].get("analysisAlias", []) + sample.get("analysisAlias", []))
                merged_samples[sample_name]["analysisAlias"] = list(combined_aliases)

        master_metadata["sample"] = list(merged_samples.values())
        master_metadata["files"] = [f for f in master_metadata["files"] if f.get("fileName") != ""]

    def convert_individual_gvf(self, assembly_path, eva_retriever, individual_gvf, json_eva):
        """Converts and unpdates metadata for a single gvf file
        :param assembly_path: path to assembly
        :param eva_retriever: eva metadata
        :param individual_gvf: gvf file
        :param json_eva: JSON file
        """
        study, _, _ = self.parse_gvf_filename(individual_gvf)
        if not study:
            logger.error(f"Could not parse study from: {individual_gvf}")
            return
        study_accession = study.split("_")[0]

        study_submission_directory = os.path.join(self.base_output_dir, "submission", study)
        os.makedirs(study_submission_directory, exist_ok=True)

        base_name = os.path.basename(individual_gvf).replace(".gvf", "")
        individual_vcf_output = os.path.join(study_submission_directory, f"{base_name}.vcf")

        if individual_gvf and individual_vcf_output and assembly_path and self.project_paths:
            convert(
                gvf_input=individual_gvf,
                vcf_output=individual_vcf_output,
                assembly=assembly_path,
                paths=self.project_paths
            )
            if eva_retriever:
                print("updating JSON with vcf info")
                eva_update_metadata_with_vcf(
                    eva_retriever=eva_retriever,
                    json_eva=json_eva,
                    vcf_output=individual_vcf_output,
                    study_accession=study_accession
                )
        else:
            logger.error("Missing parameters to convert:\n\t"
                         f"gvf {individual_gvf}\n\t"
                         f"vcf {individual_vcf_output}\n\t"
                         f"assembly {assembly_path}\n\t"
                         f"or project paths")
        return individual_vcf_output

    def _process_single_assembly(self, assembly_name, files_in_assembly, gvf_name_groups, study_accession):
        """ Process an assembly by organising input and output files, retrieving EVA metadata for the assembly and
        determining associated submitted and remapped files
        :params: assembly_name e.g. GRCh38
        :params:files_in_assembly: list of gvf files associated with that assembly
        :params: gvf_name_groups:  a dictionary where gvf file is key and (study, date, assembly) is tuple
        :return: assembly_path, eva_retriever, json_eva, remapped_files, submitted_files
        """
        first_file = files_in_assembly[0]
        study, date, assembly = gvf_name_groups[first_file]
        assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output = self.set_up_inputs_and_outputs(study,
                                                                                                              date,
                                                                                                              assembly_name)
        eva_retriever, dgva_retriever = self.retrieve_metadata(json_eva, json_dgva, study_accession, assembly_path,
                                                               assembly_report_path)
        remapped_files, submitted_files = None, None
        if files_in_assembly:
            remapped_files, submitted_files = self.check_files(files_in_assembly)
        else:
            logger.error(f"No files in this assembly: {assembly_name}")
        return assembly_path, eva_retriever, json_eva, remapped_files, submitted_files

    def _group_files_by_assembly(self, gvf_files):
        """ Groups files by the assembly. Returns dictionary based on assembly
        :params: gvf file
        :return:
            assembly_to_gvffilepaths - a dictionary where key is assembly and values is gvf files {"GRCh37": [gvf_file_paths], "GRCh38": [gvf_file_paths]})
            gvf_name_groups - a dictionary where gvf file is key and (study, date, assembly) is tuple
        """
        assembly_to_gvf_file_paths = defaultdict(list)
        gvf_name_groups = {}
        for gvf_file in gvf_files:
            # gvf_file = {estd1_Redon_et_al_2006}.{YYYY-MM-DD}.{Assembly}.{Submitted/Remapped}.gvf
            study, date, assembly = self.parse_gvf_filename(gvf_file)
            assembly_to_gvf_file_paths[assembly].append(gvf_file)
            gvf_name_groups[gvf_file] = (study, date, assembly)
        return assembly_to_gvf_file_paths, gvf_name_groups

    @staticmethod
    def parse_gvf_filename(gvf_file):
        """Extracts the assembly name from the GVF file.
        Expects GVF file format to be: {estd1_Redon_et_al_2006}.{YYYY-MM-DD}.{Assembly}.{Submitted/Remapped}.gvf
        :params: gvf_file
        :return: study e.g. estd1_Redon_et_al_2006; date e.g. YYYY-MM-DD; name of assembly e.g. GRCh37
        """
        if not gvf_file:
            return None, None, None
        # expect file_name = {estd1_Redon_et_al_2006}.{YYYY-MM-DD}.{Assembly}.{Submitted/Remapped}.gvf
        file_name = os.path.basename(gvf_file)
        if not (file_name.endswith(".Submitted.gvf") or file_name.endswith(".Remapped.gvf")):
            logger.warning(f"Warning: Skipping invalid GVF file name structure: {file_name}")
            return None, None, None
        # expect base_name = ['estd1_Redon_et_al_2006.2014-04-01.GRCh37', 'Remapped', 'gvf']
        base_name = file_name.rsplit(".", 2)
        # expect base_name[0] = {estd1_Redon_et_al_2006}.{YYYY-MM-DD}.{Assembly}
        parts = base_name[0].split(".", 2) # maxsplit 2 to handle patched assemblies i.e. assembly can be GRCh37 or GRCh37.p13
        if len(parts) < 3:
            logger.warning(f"Warning: Skipping invalid GVF file name structure: {file_name}")
            return None, None, None
        study = parts[0]
        date = parts[1]
        assembly = parts[2]
        return study, date, assembly

    def _reconfigure_assembly_metadata(self, json_eva, gvf_files, master_metadata, vcf_path):
        """Reconfigures the EVA metadata JSON with multiple files. This affects analysis, files and sample sections.
        :params: json_eva Path to EVA JSON file
        :params: gvf_files list of GVF files
        :params: master_metadata : to be updated
        :return: metadata - reconfigured for multiple assemblies
        """
        try:
            with open(json_eva, 'r') as f_in:
                metadata = json.load(f_in)
        except (FileNotFoundError, json.JSONDecodeError) as err:
            logger.error(f"Failed to update multi-analysis schema in JSON due to error: {err}")
            return None
        # find the important sections of the JSON file, leave the JSON file alone if not present.
        if not all(metadata.get(key) for key in ["analysis", "files", "sample"]):
            return
        # update analysis and file blocks with new analysis aliases
        new_analysis_aliases = self._update_analysis_and_file_blocks(
            metadata_to_add=metadata,
            files=gvf_files,
            master_metadata=master_metadata,
            vcf_path=vcf_path
        )
        # after going through each GVF file, update the analysis alias in the sample block with the new analysis aliases
        self._update_sample_block(
            new_analysis_aliases=new_analysis_aliases,
            metadata_to_add=metadata,
            master_metadata=master_metadata
        )
        return metadata

    @staticmethod
    def _update_sample_block(new_analysis_aliases, metadata_to_add, master_metadata):
        """Updates sample part of EVA JSON with new analysis aliases
        :params new_analysis_aliases: list of new names
        :params sample_list: block of EVA JSON
        """
        sample_list = metadata_to_add.get("sample", [])

        for initial_sample_block in sample_list:
            if not isinstance(initial_sample_block, dict):
                continue

            for alias in new_analysis_aliases:
                copied_sample = copy.deepcopy(initial_sample_block)
                if "analysisAlias" in copied_sample and isinstance(copied_sample["analysisAlias"], list):
                    copied_sample["analysisAlias"] = [alias]
                else:
                    copied_sample["analysisAlias"] = alias
                master_metadata["sample"].append(copied_sample)

    @staticmethod
    def _update_analysis_and_file_blocks(metadata_to_add, files, master_metadata, vcf_path):
        """Updates analysis and file parts of the EVA JSON
        :params metadata_to_add: metadata_to_add
        :params files: gvf files
        :params master_metadata: master metadata to be updated
        :return new_analysis_aliases
        """
        analysis_list = metadata_to_add.get("analysis", [])
        files_list = metadata_to_add.get("files", [])
        # get the first blocks of the EVA JSON schema and expand on those
        raw_analysis = analysis_list[0] if isinstance(analysis_list, list) and analysis_list else analysis_list
        raw_file = files_list[0] if isinstance(files_list, list) and files_list else files_list

        initial_analysis_block = raw_analysis if isinstance(raw_analysis, dict) else {}
        initial_file_block = raw_file if isinstance(raw_file, dict) else {}
        #
        # if isinstance(analysis_list, list):
        #     master_metadata["analysis"].extend(copy.deepcopy(analysis_list))
        # if isinstance(files_list, list):
        #     master_metadata["files"].extend(copy.deepcopy(files_list))

        new_analysis_aliases = []

        # change the Analysis and File blocks for each GVF file
        for index, file_path in enumerate(files, start=1):

            # rename the analysis alias string
            base_alias = initial_analysis_block.get("analysisAlias", "analysis")
            unique_alias = f"{base_alias}_file_{index}"
            new_analysis_aliases.append(unique_alias)

            # copy the analysis block and update that with new analysis alias
            analysis_block = copy.deepcopy(initial_analysis_block)
            analysis_block["analysisAlias"] = unique_alias
            master_metadata["analysis"].append(analysis_block)

            # copy the file block and update that with new analysis alias
            file_block = copy.deepcopy(initial_file_block)
            file_block["analysisAlias"] = unique_alias
            file_block["fileName"] = vcf_path
            master_metadata["files"].append(file_block)
        return new_analysis_aliases

    def retrieve_metadata(self, json_eva, json_dgva, study_accession, assembly_path, assembly_report_path):
        """Retrieves metadata and if applicable, returns list of submitted and remapped files
        :params  json_eva, json_dgva, study_accession, assembly_path, assembly_report_path
        :returns eva_retriever, dgva_retriever or eva_retriever, dgva_retriever, submitted_files, remapped_files
        """
        # fetch the metadata
        eva_retriever, dgva_retriever = gather_metadata_workflow(
            config=self.config_path,
            json_eva=json_eva,
            json_dgva=json_dgva,
            study_accession=study_accession,
            assembly=assembly_path,
            assembly_report=assembly_report_path
        )

        return eva_retriever, dgva_retriever

    @staticmethod
    def check_files(files_in_assembly):
        # for each assembly version, check if multiple submitted files
        submitted_files = [f for f in files_in_assembly if "Submitted" in os.path.basename(f)]
        remapped_files = [f for f in files_in_assembly if "Remapped" in os.path.basename(f)]
        return remapped_files, submitted_files

    def set_up_inputs_and_outputs(self, study, date, assembly_name):
        """Programmatically sets up input and output paths.
        :params assembly_name e.g. GRCh37
        :params study_accession e.g. estd1
        :return assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output: input, input, output, output, output
        """
        # get the inputs
        assembly_path = self.project_paths.assembly_paths.get(assembly_name)
        assembly_report_path = self.project_paths.assembly_report_paths.get(assembly_name)

        # programmatically generate the output files by assembly
        json_eva = os.path.join(self.base_output_dir, f"{study}.{date}.{assembly_name}.eva.json")
        json_dgva = os.path.join(self.base_output_dir, f"{study}.{date}.{assembly_name}.dgva.json")
        vcf_output = os.path.join(self.base_output_dir, f"{study}.{date}.{assembly_name}.vcf")
        return assembly_path, assembly_report_path, json_dgva, json_eva, vcf_output
