/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowPhageannotator.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Consisting of local modules
//

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { GENOMAD_DOWNLOAD                      } from '../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND                      } from '../modules/nf-core/genomad/endtoend/main'
include { CHECKV_DOWNLOADDATABASE               } from '../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND                       } from '../modules/nf-core/checkv/endtoend/main'
include { BLAST_MAKEBLASTDB                     } from '../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                          } from '../modules/nf-core/blast/blastn/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHAGEANNOTATOR {

    ch_versions = Channel.empty()

    // Import input samplesheet data using nf-validation
    Channel
        .fromSamplesheet("input")
        .multiMap { meta, fastq_1, fastq_2, fasta ->
            fastq: [ meta, [ fastq_1, fastq_2 ]]
            fasta: [ meta, fasta ]
        }
        .set { ch_input }

    //
    // MODULE: Run FastQC on reads
    //
    FASTQC (
        ch_input.fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    /*
    ----------------------------------------------------------
    Phage identification (de novo)
    ----------------------------------------------------------
    */
    //
    // MODULE: Download geNomad's database
    //
    if ( params.genomad_db ) {
        ch_genomad_db = file(params.genomad_db, checkIfExists: true)
    } else {
        ch_genomad_db = GENOMAD_DOWNLOAD ( ).genomad_db
        ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions.first())
    }

    //
    // MODULE: Identify viral sequences using geNomad
    //
    GENOMAD_ENDTOEND ( ch_input.fasta, ch_genomad_db )
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions.first())

    /*
    ----------------------------------------------------------
    Quality Filter
    ----------------------------------------------------------
    */
    //
    // MODULE: Download checkV's database
    //
    if ( params.checkv_db ) {
        ch_checkv_db = file(params.checkv_db, checkIfExists: true)
    } else {
        ch_checkv_db = CHECKV_DOWNLOADDATABASE ( ).checkv_db
        ch_versions = ch_versions.mix(CHECKV_DOWNLOADDATABASE.out.versions.first())
    }

    //
    // MODULE: Quality filter viral sequences with CheckV
    //
    CHECKV_ENDTOEND ( GENOMAD_ENDTOEND.out.virus_fasta, ch_checkv_db )
    ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions.first())

    /*
    ----------------------------------------------------------
    ANI-Clustering
    ----------------------------------------------------------
    */
    //
    // MODULE: Make BLAST database
    //
    blast_db = BLAST_MAKEBLASTDB ( CHECKV_ENDTOEND.out.viruses.map{ it[1] } ).db
    BLAST_BLASTN ( CHECKV_ENDTOEND.out.viruses, blast_db )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPhageannotator.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPhageannotator.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
