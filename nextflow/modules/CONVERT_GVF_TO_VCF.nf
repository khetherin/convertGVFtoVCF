process CONVERT_GVF_TO_VCF {
    tag {"Finding and converting GVF files for: ${study_accession}"}
    
    publishDir "${params.output_dir}", mode: 'copy' // copy to output directory

    input:
    val study_accession                 // can be study_accession string or null
    path input_dir                      // data directory
    path config_file                    // TEST.config
    path finder_script
    val credentials
    //path 'renamed_fasta*'
    output:
    val "conversion_done", emit: status_trigger
    
    script:
    """
    export REF_PATH="${params.clean_assembly_dir}"
    
    ${params.executable.convert_gvf.interpreter} ${finder_script} \\
        --search_dir ${input_dir} \\
        --log hpc.log \\
        --output "${params.output_dir}" \\
        --config ${config_file} \\
        --study_accession ${study_accession}
    """
}