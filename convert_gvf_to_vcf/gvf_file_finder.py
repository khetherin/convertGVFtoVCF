import argparse
from datetime import datetime
import hashlib
import os
#####
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from convert_gvf_to_vcf.gvf_metadata_coordinator import GvfMetadataCoordinator
from convert_gvf_to_vcf.project_paths import ProjectPaths
from convert_gvf_to_vcf.utils import get_validated_value, log_environment_packages, get_version

logger = log_cfg.get_logger(__name__)
class GvfFileFinder:
    """ The responsibility of this class is to traverse the directory to find the correct GVF file paths."""
    def __init__(self, search_dir):
        self.paths = ProjectPaths()
        self.search_dir = search_dir

    def _find_study_dirs(self, study_accession):
        """ Gathers a list of study directories. This corresponds to level 2 of the structure: /data_dir/estd3_Name_et_al_2008/gvf
        :param: study_accession
        :return: list of study directories e.g. ['estd1_Redon_et_al_2006', 'estd3_Wang_et_al_2008']
        """
        return [
            d for d in sorted(os.listdir(self.search_dir))
            if os.path.isdir(os.path.join(self.search_dir, d))
            and not d.startswith('.')
            and d.startswith(("estd", "nstd"))
            and (study_accession is None or d.startswith(study_accession))
        ]

    def _collect_gvf_files(self, gvf_dir):
        """Gathers a list of GVF files. This corresponds with level 3 of /data_dir/estd3_Name_et_al_2008/gvf
        :param: gvf_dir: full path to GVF file
        :return: list of full paths to GVF files
        """
        all_files = [f for f in os.listdir(gvf_dir) if os.path.isfile(os.path.join(gvf_dir, f))]
        non_raw_gvf = [f for f in all_files if os.path.splitext(f)[1].lower() != ".gvf"]
        if non_raw_gvf:
            logger.warning(f"Non-GVF files found in {gvf_dir}: {non_raw_gvf}")
        return [os.path.join(gvf_dir, f) for f in all_files if os.path.splitext(f)[1].lower() == ".gvf"]

    def scan(self, study_accession=None):
        """Scans a directory and returns a dictionary of the {study_accession: [gvf_files]}
        This corresponds to level 1 of /data_dir/estd3_Name_et_al_2008/gvf
        :params: study_accession if applicable
        :return: dictionary of the {study_accession: [gvf_files]}
        """
        if study_accession:
            logger.info(f"Target study accession set to: {study_accession}")
        study_dirs = self._find_study_dirs(study_accession)
        study_and_files = {}
        for study_dir in study_dirs:
            accession = study_dir.split("_")[0]
            gvf_dir = os.path.join(self.search_dir, study_dir, "gvf")
            if os.path.isdir(gvf_dir):
                study_and_files[accession] = self._collect_gvf_files(gvf_dir)
            else:
                logger.debug(f"No GVF directory found for {accession}")
        return self.deduplicate_files(study_and_files)

    def deduplicate_files(self, study_and_files):
        """ Deduplicate the GVF files found in the dictionary of lists {study_accession: [listofGVFs]}
        :params study_and_files: dictionary of list with duplicated elements
        :returns study_and_unique_files: dictionary of list with unique elements
        """
        study_and_unique_files = {}
        for study_accession, gvf_list in study_and_files.items():
            seen_files = {}
            unique_files = []
            for gvf in gvf_list:
                file_size = os.path.getsize(gvf)
                files_with_the_same_size = seen_files.get(file_size, [])
                matched_file = next((existing_file for existing_file in files_with_the_same_size if self._check_md5_match(gvf, existing_file)), None)
                if matched_file:
                    logger.error(f"Duplicate files have been found: {gvf} and {matched_file}")
                else:
                    unique_files.append(gvf)
                    seen_files[file_size] = files_with_the_same_size + [gvf]
            study_and_unique_files[study_accession] = unique_files
        return study_and_unique_files

    @staticmethod
    def get_md5(file_path):
        """Gets md5 of a file path
        :params: file_path: gvf file
        :returns: md5
        """
        hasher = hashlib.md5()
        try:
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(65536), b""):
                    hasher.update(chunk)
            return hasher.hexdigest()
        except (OSError, IOError):
            return None

    def _check_md5_match(self, file_path_1, file_path_2):
        """Returns True if both files have the exact same MD5 hash."""
        md5_path1 = self.get_md5(file_path_1)
        md5_path2 = self.get_md5(file_path_2)
        if md5_path1 and md5_path2:
            return md5_path1 == md5_path2
        return False



def main():
    VERSION = get_version()
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--search_dir", required=True, help="Directory to search i.e. FTP dir.")
    parser.add_argument("--log", required=True, help="Path to log")
    parser.add_argument("--study_accession", help="Study accession.")
    #########################################
    # for gvf file co ordination only
    parser.add_argument("--output", help="Output to gvf file coordination")
    parser.add_argument("--config", help="Config")
    #########################################
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()

    # Set up logging functionality
    if args.log:
        log_cfg.add_file_handler(args.log)
        logger.info(f"The log file is {args.log}")
        file_name = __file__
        log_environment_packages(logger, file_name, VERSION)

    finder = GvfFileFinder(search_dir= args.search_dir)
    gvf_data = finder.scan(args.study_accession)
    logger.info(f"Total number of studies found: {len(gvf_data)}")
    logger.info(f"Total number of GVF files: {sum(len(files) for files in gvf_data.values())}")
    logger.info(f"Searching for GVF files complete: {args.search_dir}")

    #########################################
    # for gvf file co ordination only
    GvfMetadataCoordinator(gvf_data, args.output, args.config).process_studies()
    #TODO: count studies stats and top dir stats
    #########################################


if __name__ == "__main__":
    main()