#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// Include modules
include { PARSE_CREDENTIALS } from './modules/PARSE_CREDENTIALS.nf'
include { CONVERT_GVF_TO_VCF } from './modules/CONVERT_GVF_TO_VCF.nf'
include { GET_ASSEMBLY_PATHS } from './modules/GET_ASSEMBLY_PATHS.nf'
include { RENAME_CONTIGS } from './modules/RENAME_CONTIGS.nf'
include { VALIDATE_SUBMISSION } from './modules/VALIDATE_SUBMISSION.nf'
include { SUBMIT_TO_EVA } from './modules/SUBMIT_TO_EVA.nf'

workflow {

    log.info """
        =======================================================
        ConvertGVFtoVCF Nextflow Pipeline Startup
        =======================================================
        Study Accession   : ${params.study_accession}
        Input Directory   : ${params.input_dir}
        Tool Config File  : ${params.tool_config}
        Output Directory  : ${params.output_dir}
        =======================================================
    """
    // Step 1: SET UP CHANNELS
    //VALUE CHANNELS - static
    input_dir_ch     = Channel.value(file(params.input_dir, type: 'dir', checkIfExists: true))
    config_file_ch   = Channel.value(file(params.tool_config, checkIfExists: true))
    credentials_ch = PARSE_CREDENTIALS(config_file_ch)
    finder_script_ch = Channel.value(file("${params.executable.convert_gvf.script_path}/gvf_file_finder.py", checkIfExists: true))
    // QUEUE CHANNELS - dynamic and trigger a new parallel task down stream
    gvf_files_ch   = Channel.fromPath("${params.input_dir}/*/gvf/*.gvf")
    // Step 2: for each GVF file prints the following to the work dir: assembly_name, fasta, report, genbank_accession
    GET_ASSEMBLY_PATHS(gvf_files_ch)

    assembly_ch = GET_ASSEMBLY_PATHS.out.map { gvf, assembly_file, fasta_file, report_file, accession_file ->
        def fasta_str = fasta_file.text.trim()
        def accession_str = accession_file.text.trim()

        def match = (fasta_str =~ /\/([^\/]+)\/[^\/]+\/[^\/]+\.[a-zA-Z0-9]+$/)
        def species_name = match.find() ? match[0][1] : "unknown_species"

        // KEY: this is the gvf and its assembly, fa, assembly report and genbank_accession
        return tuple(
            gvf,
            assembly_file.text.trim(),
            fasta_str,
            report_file.text.trim(),
            accession_str,
            species_name
        )
    }
    // Step 3 : ENSURE CONSISTENT CHROMOSOME NAMING CONVENTION FOR THE ASSEMBLY
    // RENAME_CONTIGS creates in the work dir: fasta and assembly report
    // then it copies those files (fasta and assembly report) to the output/clean_reference_sequences/species/genbank_accession
    RENAME_CONTIGS(assembly_ch)

    // Step 4 : CONVERT GVF TO VCF
    study_accession_ch = params.study_accession ? Channel.value(params.study_accession) : gvf_files_ch.map { file -> file.name.tokenize('_')[0] }.unique()

    CONVERT_GVF_TO_VCF(
        study_accession_ch,
        input_dir_ch,
        config_file_ch,
        finder_script_ch,
        credentials_ch
        //,
        //RENAME_CONTIGS.out.renamed_fasta.collect()
    )
    /*
    // Step 5: Validate submission
    study_names_ch = CONVERT_GVF_TO_VCF.out.status_trigger
        .flatMap { token ->
            def pattern = params.study_accession ? "${params.output_dir}/submission/${params.study_accession}*" : "${params.output_dir}/submission/{e,n}std[0-9]*_*"
            return file(pattern, type: 'dir').collect { it.name }
        }

    study_accessions_ch = study_names_ch.map { name -> name.split('_')[0] }

    VALIDATE_SUBMISSION(
        study_names_ch,
        study_accessions_ch,
        CONVERT_GVF_TO_VCF.out.status_trigger
    )
    // Step 6: Submit submission to EVA
    successful_logs_ch = VALIDATE_SUBMISSION.out.validation_log
        .filter { log_file ->
            log_file.readLines().any { line -> line.contains("Validation result: SUCCESS") }
        }
    empty_logs_ch = successful_logs_ch.ifEmpty("No studies passed validation successfully. Skipping SUBMIT_TO_EVA.")
    submit_inputs = successful_logs_ch.multiMap { log_file ->
        log:  log_file
        dir:  "${params.output_dir}/submission/${log_file.parent.name}"
        json: file("${params.output_dir}/submission/${log_file.parent.name}/eva_submission_*.json")
    }
    SUBMIT_TO_EVA(successful_logs_ch, PARSE_CREDENTIALS.out.credentials, submit_inputs.dir, submit_inputs.json)
    */
}