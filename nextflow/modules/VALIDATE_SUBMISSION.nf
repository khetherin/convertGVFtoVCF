process VALIDATE_SUBMISSION {
    tag "Validating submission data: ${input_study_name}"
    
    publishDir { "${params.output_dir}/validation_reports/${input_study_name}" }, mode: 'copy'

    input:
    tuple val(study_accession), val(input_study_name), val(trigger_token)

    output:
    path "**/*.json", emit: validation_manifests, optional: true
    path ".command.log", emit: validation_log

    script:
    def submit_study_dir = "${params.output_dir}/submission/${input_study_name}"
    def metadata_json = "${submit_study_dir}/eva_submission_${study_accession}.json"
    """
    source ${params.executable.eva_sub_cli.script_path}
    eva-sub-cli.py \\
        --submission_dir "${submit_study_dir}" \\
        --metadata_json "${metadata_json}" \\
        --tasks validate \\
        --validation_tasks ${params.VALIDATION_TASKS}
    """
}
