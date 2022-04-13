#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def json_schema = "$projectDir/nextflow_schema.json"

if (params.help){
  def command = "nextflow run peterk87/nf-iav-illumina --input --platfrom illumina samplesheet.csv -profile singularity/docker/conda"
  log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
  exit 0
}

if (params.validate_params) {
  NfcoreSchema.validateParameters(params, json_schema, log)
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  exit 1, "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
//custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

def modules = params.modules.clone()

//=============================================================================
// LOG PARAMS SUMMARY
//=============================================================================

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

if (params.platform == 'illumina'){
    include { ILLUMINA } from './workflows/illumina'
} else if (params.platform == 'nanopore'){
    include { NANOPORE } from './workflows/nanopore'
}

workflow AVIAN_INFLUENZA {
    if (params.platform == 'illumina'){
        ILLUMINA ()
    } else if (params.platform == 'nanopore') {
        NANOPORE ()
    }
}

workflow {
    AVIAN_INFLUENZA ()
}
/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  println """
  Pipeline execution summary
  ---------------------------
  Completed at : ${workflow.complete}
  Duration     : ${workflow.duration}
  Success      : ${c_bold}${workflow.success ? c_green : c_red}${workflow.success}${c_reset}
  Results Dir  : ${file(params.outdir)}
  Work Dir     : ${workflow.workDir}
  Exit status  : ${workflow.exitStatus}
  Error report : ${workflow.errorReport ?: '-'}
  """.stripIndent()
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
