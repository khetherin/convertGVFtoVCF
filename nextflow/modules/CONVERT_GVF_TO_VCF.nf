process CONVERT_GVF_TO_VCF {
    tag {"Finding and converting GVF files for: ${gvf_simple_name.split('_')[0]}"}
    
    publishDir "${params.output_dir}", mode: 'copy' // copy to output directory

    input:
    path input_dir                      // data directory
    path config_file                    // TEST.config
    path finder_script                  // path to gvf_file_finder.py
    val credentials                     // parsed credentials
    // e.g. homo_sapiens, GCA_000001405.1, estd1_Surname_et_al_2006, renamed_assembly(wrt GVF chromosome naming convention)
    tuple val(species), val(assembly_accession), val(gvf_simple_name), path(renamed_fasta)
    output:
    //study_accession, study_name, conversion_done
    tuple val({gvf_simple_name.split('_')[0]}), val(gvf_simple_name), val("conversion_done"), emit: status_trigger
    
    script:
    def study_accession = gvf_simple_name.split('_')[0]
    """
    export REF_PATH="${params.clean_assembly_dir}/${gvf_simple_name}"
    
    ${params.executable.convert_gvf.interpreter} ${finder_script} \\
        --search_dir ${input_dir} \\
        --log hpc.log \\
        --output "${params.output_dir}" \\
        --config ${config_file} \\
        --study_accession ${study_accession}
    """
}