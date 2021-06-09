#!/usr/bin/env nextflow
/*
 *
 ===========================================
 nf-iav-illumina: Influenza genome analysis
 ===========================================
 # Homepage / Documentation
 https://github.com/peterk87/nf-iav-illumina
 # Authors
 Peter Kruczkiewicz
 -------------------------------------------
 */

nextflow.enable.dsl = 2

def json_schema = "$projectDir/nextflow_schema.json"

if (params.help){
  def command = "nextflow run peterk87/nf-iav-illumina --input samplesheet.csv -profile singularity/docker/conda"
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
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

def modules = params.modules.clone()

//=============================================================================
// LOG PARAMS SUMMARY
//=============================================================================

def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)

include { IRMA } from './modules/local/irma'
include { CHECK_SAMPLE_SHEET } from './modules/local/check_sample_sheet'
include { SUBTYPING_REPORT } from './modules/local/subtyping_report'

include { GUNZIP as GUNZIP_FLU_FASTA } from './modules/nf-core/software/gunzip/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/software/blast/makeblastdb/main' addParams( options: modules['blastn_makeblastdb'] )
include { BLAST_BLASTN } from './modules/nf-core/software/blast/blastn/main' addParams( options: modules['blast_blastn'] )
include { CAT_FASTQ } from './modules/nf-core/software/cat/fastq/main' addParams( options: modules['cat_fastq'] )


workflow {

  GUNZIP_FLU_FASTA(ch_influenza_db_fasta)
  BLAST_MAKEBLASTDB(GUNZIP_FLU_FASTA.out.gunzip)

  Channel.fromPath( params.input, checkIfExists: true) \
    | CHECK_SAMPLE_SHEET \
    | splitCsv(header: ['sample', 'fastq1', 'fastq2', 'single_end'], sep: ',', skip: 1) \
    | map {
      def meta = [:]
      meta.id = it.sample
      meta.single_end = it.single_end.toBoolean()
      def reads = []
      def fastq1 = file(it.fastq1)
      def fastq2
      if (!fastq1.exists()) {
        exit 1, "ERROR: Please check input samplesheet. FASTQ file 1 '${fastq1}' does not exist!"
      }
      if (meta.single_end) {
        reads = [fastq1]
      } else {
        fastq2 = file(it.fastq2)
        if (!fastq2.exists()) {
          exit 1, "ERROR: Please check input samplesheet. FASTQ file 2 '${fastq2}' does not exist!"
        }
        reads = [fastq1, fastq2]
      }
      [ meta, reads ]
    } \
    | groupTuple(by: [0]) \
    | branch { meta, reads ->
      single: reads.size() == 1
        return [ meta, reads.flatten() ]
      multiple: reads.size() > 1
        return [ meta, reads.flatten() ]
    } \
    | set { ch_input }

  // Credit to nf-core/viralrecon. Source: https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/workflows/illumina.nf#L221
  // Concatenate FastQ files from same sample if required
  CAT_FASTQ(ch_input.multiple)
    .mix(ch_input.single)
    .set { ch_cat_reads }

  IRMA(ch_cat_reads)

  BLAST_BLASTN(IRMA.out.consensus, BLAST_MAKEBLASTDB.out.db)

  ch_blast = BLAST_BLASTN.out.txt.collect({ it[1] })
  SUBTYPING_REPORT(ch_influenza_metadata, ch_blast)
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
