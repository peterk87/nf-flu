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

outdir = params.outdir

def helpMessage() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_block = params.monochrome_logs ? '' : "\033[3m";
  c_ul = params.monochrome_logs ? '' : "\033[4m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_bul = c_bold + c_ul;
  log.info"""
  =${c_dim}==========================================${c_reset}
  ${c_blue+c_bold}peterk87/nf-iav-illumina${c_reset}  ~  version ${c_purple}${workflow.manifest.version}${c_reset}
  ${c_dim}===========================================${c_reset}

    ${c_ul}Git info:${c_reset} $workflow.repository - $workflow.revision [$workflow.commitId]

  ${c_bul}Usage:${c_reset}

  The typical command for running the pipeline is as follows:
  
    nextflow run peterk87/nf-iav-illumina ${c_red}--reads "${params.reads}"${c_reset} ${c_green}--outdir ${file(params.outdir)}${c_reset}
  
  ${c_yellow+c_bold+c_block}NOTE:${c_reset} ${c_yellow}Please ensure you have ${c_bul}Singularity${c_yellow} installed prior to running this workflow.${c_reset} ${c_dim}(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)${c_reset}

  ${c_bul}Mandatory Options:${c_reset}
    ${c_red}--reads${c_reset}           Input paired-end Illumina FASTQ reads; ${c_yellow}quotes required!${c_reset} (default: ${c_red}"$params.reads"${c_reset})
  ${c_bul}Other Options:${c_reset}
    ${c_green}--outdir${c_reset}          Output results directory (default: ${c_green}${file(params.outdir)}${c_reset})
    -w/--work-dir     The temporary directory where intermediate data will be saved (default: $workflow.workDir)
    -profile          Configuration profile to use [standard, slurm] (default: "$workflow.profile")
  ${c_bul}Cluster Options:${c_reset}
    --slurm_queue     Name of SLURM queue to run workflow on; use with ${c_dim}-profile slurm${c_reset}
  """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
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


// Header log info
log.info """=======================================================
${workflow.manifest.name} v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']    = workflow.manifest.name
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = file(params.outdir)
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if (workflow.profile == 'slurm') summary['SLURM Queue'] = params.slurm_queue
log.info summary.collect { k,v -> "${k.padRight(16)}: $v" }.join("\n")
log.info "========================================="


/**
 * Get number of files in a Nextflow list of files.
 */
def n_files(files) {
  return (files.getClass() == nextflow.util.BlankSeparatedList) ? files.size() : (files.size() > 0 ? 1 : 0)
}
// ==============
// PIPELINE START
// ==============
// Get Illumina read pairs stripping away characters that don't correspond to the sample name
Channel
  .fromFilePairs( params.reads, flat: true, checkIfExists: true)
  .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1], it[2] ] }
  .ifEmpty { exit 1, "No reads specified! Please specify where your reads are, e.g. '--reads \"/path/to/reads/*R{1,2}*.fastq.gz\"' (quotes around reads path required if using `*` and other characters expanded by the shell!)"}
  .into { ch_reads; ch_reads2 }

// Assembly by IRMA
process irma {
  tag "$sample_id"
  publishDir "$outdir/irma/", pattern: "*.tar.gz", mode: 'copy'
  publishDir "$outdir/irma/$sample_id", pattern: "**/READ_COUNTS.txt", mode: 'copy', saveAs: { "${sample_id}-read-counts.tsv" }
  publishDir "$outdir/consensus/$sample_id", pattern: "**/*.fa", saveAs: { filename -> file(filename).name }, mode: 'copy'

  input:
  set val(sample_id), file(reads1), file(reads2) from ch_reads
  
  output:
  set val(sample_id), file("$sample_id/amended_consensus/*.fa") optional true into ch_consensus_seq
  set val(sample_id), file("$sample_id/tables/READ_COUNTS.txt") optional true
  set val(sample_id), file("${sample_id}.tar.gz")

  script:
  """
  IRMA FLU $reads1 $reads2 $sample_id
  tar czf ${sample_id}.tar.gz $sample_id
  """
}
// Join original input sample reads channel with IRMA consensus seq channel to
// to determine which samples had no output
ch_reads2
  .join(ch_consensus_seq, remainder: true)
  .into { ch_reads_consensus; ch_reads_consensus2 }

// Generate a report of which samples produced the expected number of genome segments
process irma_report {
  tag "$sample_id"

  input:
  set val(sample_id), file(reads1), file(reads2), file(consensus_seqs) from ch_reads_consensus2

  output:
  stdout ch_consensus_report

  script:
  n_seg = n_files(consensus_seqs)
  is_complete = n_files(consensus_seqs) == 8
  """
  echo -e "${sample_id}\t${is_complete}\t$n_seg"
  """
}

/*
ch_consensus_report
  .collectFile()
  .set { ch_collected_report }
*/

// Compile basic tab-delimited table of number of genome segments generated for 
// each sample
process compile_irma_reports {
  publishDir "$outdir", mode: 'copy'

  input:
  file('report') from ch_consensus_report.collectFile()

  output:
  file('irma-consensus-report.tsv')

  """
  echo "Sample\tComplete\tSegments" > irma-consensus-report.tsv
  cat report >> irma-consensus-report.tsv
  """
}

// Download all publicly available NCBI Influenza sequences from the NCBI FTP 
// site. Genome metadata also downloaded for report generation.
process download_influenza_db {
  output:
  file('influenza.fna') into ch_influenza_db_fasta
  file('genomeset.dat') into ch_influenza_db_metadata

  """
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz" 2> /dev/null
  gunzip influenza.fna.gz
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz" 2> /dev/null
  gunzip genomeset.dat.gz
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

// Nucleotide BLAST all segment consensus sequences for each sample against the 
// NCBI Influenza sequences.
process blastn_irma_consensus_seqs {
  tag "$sample_id"

  input:
  file(blastn_db) from ch_influenza_blastn_db
  set val(sample_id), file(reads1), file(reads2), file(consensus_seqs) from ch_reads_consensus

  output:
  set val(sample_id), file(blast_out) into ch_blast

  when:
  consensus_seqs.size() > 0

  script:
  db_name = blastn_db[0].baseName
  blast_out = "${sample_id}-vs-ncbi-influenza.tsv"
  """
  cat $consensus_seqs > segs.fa
  blastn \\
   -db $db_name \\
   -query segs.fa \\
   -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle" \\
   -out $blast_out \\
   -num_threads ${task.cpus} \\
   -num_alignments 1000000 \\
   -evalue 1e-6
  """
}

ch_blast
  .collect { it[1] }
  .set { ch_all_blastn_results }

process subtyping_report {
  publishDir "$outdir/", mode: 'copy'
  tag "$mem_reqs GB"
  memory {
   "${mem_reqs} GB" 
 }

  input:
  file('genomeset.dat') from ch_influenza_db_metadata
  file(blastn_results) from ch_all_blastn_results

  output:
  file('subtyping_report.xlsx')

  script:
  mem_reqs = Math.ceil(2 * (blastn_results.collect { it.size() }.sum()) / (1024**3)) + 2
  """
  parse_influenza_blast_results.py \\
   --threads ${task.cpus} \\
   --flu-metadata genomeset.dat \\
   --excel-report subtyping_report.xlsx \\
   $blastn_results
  """
}


/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_bold = params.monochrome_logs ? '' : "\033[1m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_black = params.monochrome_logs ? '' : "\033[0;30m";
  c_red = params.monochrome_logs ? '' : "\033[0;31m";
  c_green = params.monochrome_logs ? '' : "\033[0;32m";
  c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_purple = params.monochrome_logs ? '' : "\033[0;35m";
  c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
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
