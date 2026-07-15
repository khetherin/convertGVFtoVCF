import argparse
import json
import os.path
import shutil

from convert_gvf_to_vcf.metadata_retrievers.dgva_metadata import DGVAMetadataRetriever
from convert_gvf_to_vcf.metadata_retrievers.eva_metadata import EVAMetadataRetriever
from ebi_eva_common_pyutils.logger import logging_config as log_cfg
logger = log_cfg.get_logger(__name__)

import json
import os
import shutil


def eva_add_file_metadata(retriever, json_output, vcf_output, study_accession):
    files_analysis_id_list = retriever._fetch_analysis_ids(study_accession)
    files_analysis_alias = retriever._fetch_analysis_alias(study_accession, files_analysis_id_list)

    files_file_name = retriever._get_file_name(vcf_output)
    files_file_size = retriever._get_file_size(vcf_output)
    files_file_md5 = retriever._get_file_md5(vcf_output)

    if not os.path.exists(json_output):
        logger.error(f"Cannot find json output: {json_output}")
    # moving the file to a preconversion file to prevent confusion
    base_path, ext = os.path.splitext(json_output)
    preconversion_json_path = f"{base_path}_preconverted{ext}"

    try:
        shutil.copy(json_output, preconversion_json_path)
        with open(json_output, 'r') as f_in:
            metadata = json.load(f_in)
        if "files" not in metadata or not isinstance(metadata["files"], list):
            metadata["files"] = []
        file_found = False
        for file_object in metadata["files"]:
            if file_object.get("fileName") == files_file_name:
                file_object["analysisAlias"] = files_analysis_alias
                file_object["fileSize"] = files_file_size
                file_object["md5"] = files_file_md5
                file_found = True
                break  # Exit loop immediately once match is handled
        # If a new file, append it
        if not file_found:
            metadata["files"].append({
                "analysisAlias": files_analysis_alias,
                "fileName": files_file_name,
                "fileSize": files_file_size,
                "md5": files_file_md5
            })
        with open(json_output, 'w') as f_out:
            json.dump(metadata, f_out, indent=4)
    except (IOError, json.JSONDecodeError) as e:
        logger.error(f"Failed to process or rewrite metadata for {json_output}: {e}")


def gather_metadata_workflow(config, json_eva, json_dgva, study_accession, assembly, assembly_report):
    eva_retriever = None
    dgva_retriever = None
    if json_eva:
        eva_retriever = EVAMetadataRetriever(
            path_to_config_yaml=config
        )
        eva_retriever.create_json_eva(json_eva, study_accession, assembly, assembly_report)

    if json_dgva:
        dgva_retriever = DGVAMetadataRetriever(
            path_to_config_yaml=config
        )
        dgva_retriever.create_json_dgva(json_dgva, study_accession)
    return eva_retriever, dgva_retriever

def eva_update_metadata_with_vcf(eva_retriever, json_eva, vcf_output, study_accession):
    if eva_retriever and vcf_output and json_eva:
        eva_add_file_metadata(eva_retriever, json_eva, vcf_output, study_accession)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config")
    parser.add_argument("--json_output_eva")
    parser.add_argument("--json_output_dgva")
    parser.add_argument("--study_accession")
    parser.add_argument("--vcf_output")
    parser.add_argument("--assembly")
    parser.add_argument("--assembly_report")
    args = parser.parse_args()

    eva_retriever, dgva_retriever = gather_metadata_workflow(
        config=args.config,
        json_eva=args.json_output_eva,
        json_dgva=args.json_output_dgva,
        study_accession=args.study_accession,
        assembly=args.assembly,
        assembly_report=args.assembly_report
    )
    eva_update_metadata_with_vcf(eva_retriever, json_eva=args.json_output_eva, vcf_output=args.vcf_output, study_accession=args.study_accession)


if __name__ == "__main__":
    main()