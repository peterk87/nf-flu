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

// Default parameters
params.reads_dir = "$baseDir/"
params.reads_pattern = "*_S*_L???_R{1,2}_001.fastq.gz"
params.outdir = "$baseDir/results"
params.help = false

outdir = params.outdir

log.info "log=$log, ${log.getClass()}"
def helpMessage() {
    log.info"""
    ===================================
    peterk87/nf-iav-illumina  ~  version ${workflow.manifest.version}
    ===================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run peterk87/nf-iav-illumina --reads_dir /path/to/illumina/reads/directory/ -profile standard
    Mandatory arguments:
      --reads_dir                   Path to Illumina FASTQ file directory
    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      --reads_pattern               Glob pattern to match reads (default for raw Illumina Gzipped FASTQ ("*_S*_L???_R{1,2}_001.fastq.gz"))
      --slurm_queue                 Name of SLURM queue to run workflow on (default "")
      -profile                      Configuration profile to use. [standard, slurm] (default 'standard')
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if (workflow.profile == 'slurm' && params.slurm_queue == "") {
  log.error "You must specify a valid SLURM queue (e.g. '--slurm_queue <queue name>' (see `\$ sinfo` output for available queues)) to run this workflow with the 'slurm' profile!"
  exit 1
}

log.info """
==============================
Pipeline Initialization Info
==============================
Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
Params  : $params
Nextflow version: $workflow.nextflow.version
------------------------------
"""

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
  .fromFilePairs( params.reads_dir + params.reads_pattern, flat: true)
  .map { [ it[0].replaceAll(/_S\d{1,2}_L001/, ""), it[1], it[2] ] }
  .into { samples_ch; samples_ch2 }

// Assembly by IRMA
process irma {
  tag "$sample_id"
  publishDir "$outdir/irma/", pattern: "*.tar.gz", mode: 'copy'
  publishDir "$outdir/irma/$sample_id", pattern: "**/READ_COUNTS.txt", mode: 'copy', saveAs: { "${sample_id}-read-counts.tsv" }
  publishDir "$outdir/consensus/$sample_id", pattern: "**/*.fa", saveAs: { filename -> file(filename).name }, mode: 'copy'

  input:
  set val(sample_id), file(reads1), file(reads2) from samples_ch
  
  output:
  set val(sample_id), file("$sample_id/amended_consensus/*.fa") optional true into consensus_seq_ch
  set val(sample_id), file("$sample_id/tables/READ_COUNTS.txt") optional true into irma_read_counts_ch
  set val(sample_id), file("${sample_id}.tar.gz") into irma_tgz_ch

  script:
  """
  IRMA FLU $reads1 $reads2 $sample_id
  tar czf ${sample_id}.tar.gz $sample_id
  """
}
// Join original input sample reads channel with IRMA consensus seq channel to
// to determine which samples had no output
samples_ch2
  .join(consensus_seq_ch, remainder: true)
  .into { reads_consensus_ch; reads_consensus_ch2 }

// Generate a report of which samples produced the expected number of genome segments
process irma_report {
  tag "$sample_id"

  input:
  set val(sample_id), file(reads1), file(reads2), file(consensus_seqs) from reads_consensus_ch2

  output:
  stdout consensus_report_ch

  script:
  n_seg = n_files(consensus_seqs)
  is_complete = n_files(consensus_seqs) == 8
  """
  echo -e "${sample_id}\t${is_complete}\t$n_seg"
  """
}

consensus_report_ch
  .collectFile()
  .set { collected_report_ch }

// Compile basic tab-delimited table of number of genome segments generated for 
// each sample
process compile_irma_reports {
  publishDir "$outdir", mode: 'copy'

  input:
  file('report') from collected_report_ch

  output:
  file('irma-consensus-report.tsv') into compiled_consensus_report_ch

  """
  echo "Sample\tComplete\tSegments" > irma-consensus-report.tsv
  cat report >> irma-consensus-report.tsv
  """
}

// Download all publicly available NCBI Influenza sequences from the NCBI FTP 
// site. Genome metadata also downloaded for report generation.
process download_influenza_db {
  output:
  file('influenza.fna') into influenza_db_fasta
  file('genomeset.dat') into influenza_db_metadata

  """
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz" 2> /dev/null
  gunzip influenza.fna.gz
  wget "ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz" 2> /dev/null
  gunzip genomeset.dat
  """
}

// Create a nucleotide BLAST database from the Influenza sequences downloaded 
// from NCBI.
process blast_db {
  input:
  file('influenza.fna') from influenza_db_fasta

  output:
  file("*.{nhr,nin,nsq}") into influenza_blastn_db

  """
  makeblastdb -dbtype nucl -in influenza.fna
  """
}

// Nucleotide BLAST all segment consensus sequences for each sample against the 
// NCBI Influenza sequences.
process blastn_irma_consensus_seqs {
  tag "$sample_id"

  input:
  file(blastn_db) from influenza_blastn_db
  set val(sample_id), file(reads1), file(reads2), file(consensus_seqs) from reads_consensus_ch

  output:
  set val(sample_id), file(blast_out) into blastn_ch

  when:
  consensus_seqs.size() > 0

  script:
  db_name = blastn_db[0].baseName
  blast_out = "${sample_id}-vs-ncbi-influenza.tsv"
  """
  cat $consensus_seqs > segs.fa
  blastn -db $db_name -query segs.fa -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle" -out $blast_out -num_threads ${task.cpus} -num_alignments 1000000 -evalue 1e-6
  """
}

blastn_ch
  .collect { it[1] }
  .set { all_blastn_results }

process subtyping_report {
  publishDir "$outdir/", mode: 'copy'

  input:
  file('genomeset.dat') from influenza_db_metadata
  file(blastn_results) from all_blastn_results

  output:
  file('subtyping_report.xlsx') into subtyping_report_xlsx

  script:
  script_path = "$baseDir/scripts/parse_influenza_blast_results.py"

  """
  python $script_path --threads ${task.cpus} --flu-metadata genomeset.dat --excel-report subtyping_report.xlsx $blastn_results
  """
}


/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Results Dir  : ${file(params.outdir)}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """.stripIndent()
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
