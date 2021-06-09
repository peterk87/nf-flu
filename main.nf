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



include { IRMA } from './modules/local/irma'

include { GUNZIP as GUNZIP_FLU_FASTA } from './modules/nf-core/software/gunzip/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/software/blast/makeblastdb/main' addParams( options: modules['blastn_makeblastdb'] )
include { BLAST_BLASTN } from './modules/nf-core/software/blast/blastn/main' addParams( options: modules['blast_blastn'] )


/**
 * Get number of files in a Nextflow list of files.
def n_files(files) {
  return (files.getClass() == nextflow.util.BlankSeparatedList) ? files.size() : (files.size() > 0 ? 1 : 0)
}

// Join original input sample reads channel with IRMA consensus seq channel to
// to determine which samples had no output
ch_reads2
  .join(ch_consensus_seq, remainder: true)
  .into { ch_reads_consensus; ch_reads_consensus2 }

// Generate a report of which samples produced the expected number of genome segments
process irma_report {
  tag "$sample"

  input:
  tuple val(sample), path(reads), path(consensus_seqs)

  output:
  stdout emit: stdout

  script:
  n_seg = n_files(consensus_seqs)
  is_complete = n_files(consensus_seqs) == 8
  """
  echo -e "${sample}\t${is_complete}\t$n_seg"
  """
}

/*
ch_consensus_report
  .collectFile()
  .set { ch_collected_report }


// Compile basic tab-delimited table of number of genome segments generated for 
// each sample
process compile_irma_reports {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple('report') // .collectFile()

  output:
  file('irma-consensus-report.tsv')

  """
  echo "Sample\tComplete\tSegments" > irma-consensus-report.tsv
  cat report >> irma-consensus-report.tsv
  """
}

// Create a nucleotide BLAST database from the Influenza sequences downloaded 
// from NCBI.
process blast_db {
  input:
  file('influenza.fna') from ch_influenza_db_fasta

  output:
  file("*.{nhr,nin,nsq}") into ch_influenza_blastn_db

  """
  makeblastdb -dbtype nucl -in influenza.fna
  """
}


ch_blast
  .collect { it[1] }
  .set { ch_all_blastn_results }
*/
process SUBTYPING_REPORT {
  publishDir "${params.outdir}/", mode: params.publish_dir_mode
  memory { 
    // Dynamically determine how much memory is required for this task based on 
    // overall size of tabular blastn inputs. For a single input, allocate 2GB
    if (blastn_results instanceof nextflow.processor.TaskPath) {
      // single input file allocate 2GB
      "2 GB"
    } else if (blastn_results instanceof nextflow.util.BlankSeparatedList) {
      // multiple input files
      // mem reqs = half of sum of GB file sizes plus 2GB wiggle room
      mem_reqs = Math.ceil(0.5 * (blastn_results.collect { it.size() }.sum()) / (1024**3)) + 2
      "${mem_reqs} GB" 
    } else {
      // not TaskPath or BlankSeparatedList, then default 2GB for memory
      "2 GB"
    }
  }

  input:
  path(genomeset)
  path(blastn_results)

  output:
  path('iav-subtyping-report.xlsx'), emit: report
  path('parse_influenza_blast_results.log'), emit: log

  script:
  """
  parse_influenza_blast_results.py \\
   --threads ${task.cpus} \\
   --flu-metadata $genomeset \\
   --excel-report iav-subtyping-report.xlsx \\
   --pident-threshold $params.pident_threshold \\
   $blastn_results
  ln -s .command.log parse_influenza_blast_results.log
  """
}


process CHECK_SAMPLE_SHEET {
  publishDir "${params.tracedir}/",
             mode: params.publish_dir_mode
  // using shiptv container since it has pandas, rich, typer installed
  conda (params.enable_conda ? 'bioconda::shiptv=0.4.0' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  path(samplesheet)

  output:
  path('samplesheet.fixed.csv')

  script:
  """
  check_sample_sheet.py $samplesheet samplesheet.fixed.csv
  """
}


def check_sample_sheet(LinkedHashMap sample_sheet) {
  // Check that each entry from a sample sheet
  reads_file = file(sample_sheet.reads, checkIfExists: true)
  reads = sample_sheet.reads ? reads_file : null
  if (reads == null) {
    exit 1, "The Nanopore reads FASTQ file or directory specified for ${sample_sheet.sample} does not exist! Please check the sample sheet '$params.sample_sheet'"
  }
  if (reads_file.isDirectory()) {
    fqs = []
    reads_file.eachFile {
      fname = it.getName()
      if (fname ==~ /.*\.fastq(\.gz)?$/) {
        fqs << it
      }
    }
    if (fqs.size() == 0) {
      exit 1, "Sample '${sample_sheet.sample}' input specified as directory '$reads_file' with no FASTQ files! Please check the sample sheet '$params.sample_sheet'\n${fqs}\n${reads_file.listFiles()}"
    }
    reads = fqs
  }
  return [ sample_sheet.sample, reads ]
}



workflow {
  // Get Illumina read pairs stripping away characters that don't correspond to the sample name
  /*
  ch_reads = Channel.fromFilePairs( params.reads, checkIfExists: true)
    .map { [ [id: it[0].replaceAll(/_S\d{1,2}_L001/, "")], it[1] ] }
    .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}*.fastq.gz\"' (quotes around reads path required if using `*` and other characters expanded by the shell!)"}
  */
  ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
  ch_influenza_metadata = file(params.ncbi_influenza_metadata)
  GUNZIP_FLU_FASTA(ch_influenza_db_fasta)
  BLAST_MAKEBLASTDB(GUNZIP_FLU_FASTA.out.gunzip)

  Channel.fromPath( params.input, checkIfExists: true) \
    | CHECK_SAMPLE_SHEET \
    | splitCsv(header: ['sample', 'fastq1', 'fastq2', 'single_end'], sep: ',', skip: 1) \
    | map { 
      [ [id: it.sample, single_end: it.single_end], [it.fastq1, it.fastq2] ]
    } \
    | set { ch_input }

  IRMA(ch_input)

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
