#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//=============================================================================
// NCBI Influenza DB reference data
//=============================================================================

ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)

//=============================================================================
// NCBI VADR Influenza virus model
//=============================================================================

ch_vadr_model_targz = file(params.vadr_model_targz)

//=============================================================================
// MODULES
//=============================================================================

include { SUBTYPING_REPORT } from '../modules/local/subtyping_report'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI } from '../modules/local/blast_makeblastdb'
include { BLAST_BLASTN } from '../modules/local/blastn'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_FASTA; ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_CSV } from '../modules/local/zstd_decompress'
include { SETUP_FLU_VADR_MODEL } from '../modules/local/vadr'
include { VADR; VADR_SUMMARIZE_ISSUES } from '../modules/local/vadr'
include { VADR2BLASTN } from '../modules/local/vadr'
include { PRE_TABLE2ASN } from '../modules/local/table2asn'
include { POST_TABLE2ASN } from '../modules/local/table2asn'
include { TABLE2ASN } from '../modules/local/table2asn'
include { MQC_VERSIONS_TABLE } from '../modules/local/mqc_versions_table'
include { FLUMUT; VADR2FLUMUT } from '../modules/local/flumut'
include { GENOFLU } from '../modules/local/genoflu'
include { CLEAVAGE_SITE } from '../modules/local/cleavage_site'

//=============================================================================
// Workflow Params Setup
//=============================================================================

def summary_params = NfcoreSchema.params_summary_map(workflow, params, "$projectDir/nextflow_schema.json")

//=============================================================================
// WORKFLOW
//=============================================================================
workflow ASSEMBLIES {
  ch_versions = Channel.empty()

  ch_input_fasta = Channel.fromPath("${params.input}/*.fa*").map {
    def baseName = it.baseName
    // remove certain suffixes from the base filename
    baseName = baseName.replaceAll(/\.consensus$/, '')
    [baseName, it] 
  }

  // Decompress reference data
  ZSTD_DECOMPRESS_FASTA(ch_influenza_db_fasta, "influenza.fasta")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_FASTA.out.versions)
  ZSTD_DECOMPRESS_CSV(ch_influenza_metadata, "influenza.csv")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_CSV.out.versions)
  BLAST_MAKEBLASTDB_NCBI(ZSTD_DECOMPRESS_FASTA.out.file)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_NCBI.out.versions)
  SETUP_FLU_VADR_MODEL(ch_vadr_model_targz, params.custom_flu_minfo)


  VADR(ch_input_fasta, SETUP_FLU_VADR_MODEL.out)
  ch_versions = ch_versions.mix(VADR.out.versions)
  VADR.out.feature_table
    .combine(VADR.out.pass_fasta, by: 0)
    .set { ch_pre_table2asn }
  VADR_SUMMARIZE_ISSUES(VADR.out.vadr_outdir.map { [it[1]] }.collect())
  PRE_TABLE2ASN(ch_pre_table2asn)
  ch_versions = ch_versions.mix(PRE_TABLE2ASN.out.versions)
  TABLE2ASN(PRE_TABLE2ASN.out.table2asn_input)
  ch_versions = ch_versions.mix(TABLE2ASN.out.versions)
  POST_TABLE2ASN(TABLE2ASN.out.genbank)
  ch_versions = ch_versions.mix(POST_TABLE2ASN.out.versions)


  VADR2BLASTN(POST_TABLE2ASN.out.genbank.join(VADR.out.vadr_outdir))
  BLAST_BLASTN(
    VADR2BLASTN.out.fasta.map { [[id: it[0]], it[1]] },
    BLAST_MAKEBLASTDB_NCBI.out.db
  )
  ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first().ifEmpty(null))

  ch_blast = BLAST_BLASTN.out.txt.collect({ it[1] })
  ch_vadr_outdir = VADR.out.vadr_outdir.map { [it[1]] }.collect()
  SUBTYPING_REPORT(
    ZSTD_DECOMPRESS_CSV.out.file,
    ch_blast,
    ch_vadr_outdir,
    []
  )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT.out.versions)

  GENOFLU(ch_input_fasta)
  ch_versions = ch_versions.mix(GENOFLU.out.versions)
  
  CLEAVAGE_SITE(POST_TABLE2ASN.out.cds_aa_fasta)
  ch_versions = ch_versions.mix(CLEAVAGE_SITE.out.versions)

  if (!params.skip_flumut) {
    VADR2FLUMUT(POST_TABLE2ASN.out.genbank.collect({ it[1] }))
    FLUMUT(VADR2FLUMUT.out.fasta)
    ch_versions = ch_versions.mix(FLUMUT.out.versions)
  }

}
