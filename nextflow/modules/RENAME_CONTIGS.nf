process RENAME_CONTIGS {
    tag "Renaming contigs for ${assembly_accession}"
    
    publishDir "${params.clean_assembly_dir}/${species}/${assembly_accession}/${gvf_file.simpleName}", mode: 'symlink'
    input:
    tuple path(gvf_file), val(assembly_name), val(assembly_fasta), val(assembly_report), val(assembly_accession), val(species)

    output:
    tuple val(species), val(assembly_accession), val(gvf_file.simpleName), path("${assembly_accession}.fa"), emit: ready_to_convert
    path "${assembly_accession}_assembly_report.txt", emit: safe_report
    script:
    def output_fasta = "${assembly_accession}.fa"

    """
    export PYTHONPATH="${params.executable.eva_submission.script_path}"
    export REF_PATH="${params.REF_PATH}"
    ${params.executable.eva_submission.interpreter} -m eva_submission.steps.rename_contigs_from_insdc_in_assembly \\
        --get_contig_from_vcf data \\
        --assembly_accession "${assembly_accession}" \\
        --custom_fasta "${output_fasta}" \\
        --assembly_fasta "${assembly_fasta}" \\
        --assembly_report "${assembly_report}" \\
        --vcf_files "${gvf_file}"
    cp "${assembly_report}" "${assembly_accession}_assembly_report.txt"
    """
}
