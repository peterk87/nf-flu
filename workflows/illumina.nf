#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//=============================================================================
// NCBI Influenza DB reference data
//=============================================================================

ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)

//=============================================================================
// MODULES
//=============================================================================

include { IRMA } from '../modules/local/irma'
include { CHECK_SAMPLE_SHEET } from '../modules/local/check_sample_sheet'
include { SUBTYPING_REPORT } from '../modules/local/subtyping_report'
include { BLAST_MAKEBLASTDB } from '../modules/local/blast_makeblastdb'
include { BLAST_BLASTN } from '../modules/local/blastn'
include { CAT_ILLUMINA_FASTQ } from '../modules/local/cat_illumina_fastq'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_FASTA; ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_CSV } from '../modules/local/zstd_decompress'

include { CUSTOM_DUMPSOFTWAREVERSIONS  as SOFTWARE_VERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

//=============================================================================
// Workflow Params Setup
//=============================================================================

def irma_module = 'FLU-utr'
if (params.irma_module) {
    irma_module = params.irma_module
}

//=============================================================================
// WORKFLOW
//=============================================================================

workflow ILLUMINA {
  ch_versions = Channel.empty()
  // Decompress reference data
  ZSTD_DECOMPRESS_FASTA(ch_influenza_db_fasta, "influenza.fasta")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_FASTA.out.versions)
  ZSTD_DECOMPRESS_CSV(ch_influenza_metadata, "influenza.csv")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_CSV.out.versions)
  BLAST_MAKEBLASTDB(ZSTD_DECOMPRESS_FASTA.out.file)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

  CHECK_SAMPLE_SHEET(Channel.fromPath( params.input, checkIfExists: true))
    .splitCsv(header: ['sample', 'fastq1', 'fastq2', 'single_end'], sep: ',', skip: 1)
    .map {
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
    } 
    .groupTuple(by: [0]) \
    .map { meta, reads ->
      return [ meta, reads.flatten() ]
    }
    .set { ch_input }

  // Credit to nf-core/viralrecon. Source: https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/workflows/illumina.nf#L221
  // Concatenate FastQ files from same sample if required
  CAT_ILLUMINA_FASTQ(ch_input)
  ch_versions = ch_versions.mix(CAT_ILLUMINA_FASTQ.out.versions.first().ifEmpty(null))

  IRMA(CAT_ILLUMINA_FASTQ.out.reads, irma_module)
  ch_versions = ch_versions.mix(IRMA.out.versions.first().ifEmpty(null))

  BLAST_BLASTN(IRMA.out.consensus, BLAST_MAKEBLASTDB.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first().ifEmpty(null))

  ch_blast = BLAST_BLASTN.out.txt.collect({ it[1] })
  SUBTYPING_REPORT(
    ZSTD_DECOMPRESS_CSV.out.file,
    ch_blast,
    CHECK_SAMPLE_SHEET.out
  )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT.out.versions)

  SOFTWARE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}
