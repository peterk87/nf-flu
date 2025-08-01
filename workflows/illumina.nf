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

include { IRMA                                                    } from '../modules/local/irma'
include { CHECK_SAMPLE_SHEET                                      } from '../modules/local/check_sample_sheet'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_IRMA_CONSENSUS     } from '../modules/local/subtyping_report'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_BCF_CONSENSUS      } from '../modules/local/subtyping_report'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI             } from '../modules/local/blast_makeblastdb'
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                       } from '../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS                  } from '../modules/local/blastn'
include { CAT_ILLUMINA_FASTQ                                      } from '../modules/local/cat_illumina_fastq'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_FASTA                } from '../modules/local/zstd_decompress'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_CSV                  } from '../modules/local/zstd_decompress'
include { MULTIQC                                                 } from '../modules/local/multiqc'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV            } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV            } from '../modules/local/multiqc_tsv_from_list'
include { MINIMAP2                                                } from '../modules/local/minimap2'
include { MOSDEPTH_GENOME                                         } from '../modules/local/mosdepth'
include { BCFTOOLS_STATS                                          } from '../modules/local/bcftools'
include { SEQTK_SEQ                                               } from '../modules/local/seqtk_seq'
include { PULL_TOP_REF_ID                                         } from '../modules/local/pull_top_ref_id'
include { BCF_FILTER as BCF_FILTER_FREEBAYES                      } from '../modules/local/bcftools'
include { BCF_CONSENSUS                                           } from '../modules/local/bcftools'
include { COVERAGE_PLOT                                           } from '../modules/local/coverage_plot'
include { CAT_CONSENSUS                                           } from '../modules/local/misc'
include { FREEBAYES                                               } from '../modules/local/freebayes'
include { SETUP_FLU_VADR_MODEL                                    } from '../modules/local/vadr'
include { VADR as VADR_IRMA                                       } from '../modules/local/vadr'
include { VADR_SUMMARIZE_ISSUES as VADR_SUMMARIZE_ISSUES_IRMA     } from '../modules/local/vadr'
include { VADR as VADR_BCFTOOLS                                   } from '../modules/local/vadr'
include { VADR_SUMMARIZE_ISSUES as VADR_SUMMARIZE_ISSUES_BCFTOOLS } from '../modules/local/vadr'
include { PRE_TABLE2ASN as PRE_TABLE2ASN_IRMA                     } from '../modules/local/table2asn'
include { TABLE2ASN as TABLE2ASN_IRMA                             } from '../modules/local/table2asn'
include { POST_TABLE2ASN as POST_TABLE2ASN_IRMA                   } from '../modules/local/table2asn'
include { PRE_TABLE2ASN as PRE_TABLE2ASN_BCFTOOLS                 } from '../modules/local/table2asn'
include { TABLE2ASN as TABLE2ASN_BCFTOOLS                         } from '../modules/local/table2asn'
include { POST_TABLE2ASN as POST_TABLE2ASN_BCFTOOLS               } from '../modules/local/table2asn'
include { MQC_VERSIONS_TABLE                                      } from '../modules/local/mqc_versions_table'
include { FLUMUT; PREP_FLUMUT_FASTA                               } from '../modules/local/flumut'
include { GENOFLU                                                 } from '../modules/local/genoflu'
include { CLEAVAGE_SITE                                           } from '../modules/local/cleavage_site'
include { GENIN2                                                  } from '../modules/local/genin2'
// SUBWORKFLOWS
include { NEXTCLADE } from '../subworkflows/nextclade'

//=============================================================================
// Workflow Params Setup
//=============================================================================

def irma_module = 'FLU-utr'
if (params.irma_module) {
    irma_module = params.irma_module
}

def pass_sample_reads = [:]
def fail_sample_reads = [:]
def summary_params = NfcoreSchema.params_summary_map(workflow, params, "$projectDir/nextflow_schema.json")

//=============================================================================
// WORKFLOW
//=============================================================================

workflow ILLUMINA {
  ch_versions = Channel.empty()

  // Sample Sheet Check
  ch_input = CHECK_SAMPLE_SHEET(Channel.fromPath(params.input, checkIfExists: true))

    ch_input.splitCsv(header: ['sample', 'fastq1', 'fastq2', 'single_end'], sep: ',', skip: 1)
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
    .groupTuple(by: [0])
    .map { meta, reads ->
      return [ meta, reads.flatten() ]
    }
    .set { ch_input_sorted }

  // Read Count Check
  ch_input_sorted
    .map { meta, reads ->
      def count = reads.collect { it.countFastq() }.sum()
      return [ meta, reads, count ]
    }
    .branch { meta, reads, count ->
      pass: count >= params.min_sample_reads
        pass_sample_reads[meta.id] = count
        return [ "${meta.id}\t$count" ]
      fail: count < params.min_sample_reads
        fail_sample_reads[meta.id] = count
        return [ "${meta.id}\t$count" ]
    }
    .set { ch_pass_fail_read_count }

  // Report samples which have reads count < min_sample_reads
  READ_COUNT_FAIL_TSV(
    ch_pass_fail_read_count.fail.collect(),
    ['Sample', 'Read count'],
    'fail_read_count_samples'
  )
  // Report samples which have reads count >= min_sample_reads
  READ_COUNT_PASS_TSV(
    ch_pass_fail_read_count.pass.collect(),
    ['Sample', 'Read count'],
    'pass_read_count_samples'
  )
  // Keep samples which have reads count >= min_sample_reads for downstream analysis
  // Re-arrange channels to have meta map of information for sample
  // ch_input_sorted
  //   .filter { it[2] >= params.min_sample_reads }
  //   .map { meta, reads, count -> [ meta, reads ] }
  //   .set { ch_reads }

  // Decompress reference data
  ZSTD_DECOMPRESS_FASTA(ch_influenza_db_fasta, "influenza.fasta")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_FASTA.out.versions)
  ZSTD_DECOMPRESS_CSV(ch_influenza_metadata, "influenza.csv")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_CSV.out.versions)
  BLAST_MAKEBLASTDB_NCBI(ZSTD_DECOMPRESS_FASTA.out.file)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_NCBI.out.versions)
  SETUP_FLU_VADR_MODEL(ch_vadr_model_targz, params.custom_flu_minfo)

  // Use ch_input_sorted for CAT_ILLUMINA_FASTQ to ensure IRMA triggers
  CAT_ILLUMINA_FASTQ(ch_input_sorted)

  // IRMA processing
  IRMA(CAT_ILLUMINA_FASTQ.out.reads, irma_module)
  ch_versions = ch_versions.mix(IRMA.out.versions.first().ifEmpty(null))

  // VADR application on IRMA concsensus
  IRMA.out.consensus
    .map { [it[0].id, it[1]] }
    .set { ch_irma_consensus }
  VADR_IRMA(ch_irma_consensus, SETUP_FLU_VADR_MODEL.out)
  ch_versions = ch_versions.mix(VADR_IRMA.out.versions)
  VADR_IRMA.out.feature_table
    .combine(VADR_IRMA.out.pass_fasta, by: 0)
    .set { ch_pre_table2asn }
  VADR_SUMMARIZE_ISSUES_IRMA(VADR_IRMA.out.vadr_outdir.map { [it[1]] }.collect())
  PRE_TABLE2ASN_IRMA(ch_pre_table2asn)
  ch_versions = ch_versions.mix(PRE_TABLE2ASN_IRMA.out.versions)
  TABLE2ASN_IRMA(PRE_TABLE2ASN_IRMA.out.table2asn_input)
  ch_versions = ch_versions.mix(TABLE2ASN_IRMA.out.versions)
  POST_TABLE2ASN_IRMA(TABLE2ASN_IRMA.out.genbank)
  ch_versions = ch_versions.mix(POST_TABLE2ASN_IRMA.out.versions)

  BLAST_BLASTN_IRMA(IRMA.out.consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN_IRMA.out.versions.first().ifEmpty(null))

  ch_blast = BLAST_BLASTN_IRMA.out.txt.collect({ it[1] })
  ch_vadr_outdir_irma = VADR_IRMA.out.vadr_outdir.map { [it[1]] }.collect()
  SUBTYPING_REPORT_IRMA_CONSENSUS(
    ZSTD_DECOMPRESS_CSV.out.file,
    ch_blast,
    ch_vadr_outdir_irma,
    CHECK_SAMPLE_SHEET.out
  )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT_IRMA_CONSENSUS.out.versions)

  // Prepare top ncbi accession id for each segment of each sample (id which has top bitscore)
  PULL_TOP_REF_ID(BLAST_BLASTN_IRMA.out.txt, ZSTD_DECOMPRESS_CSV.out.file)
  ch_versions = ch_versions.mix(PULL_TOP_REF_ID.out.versions)

  PULL_TOP_REF_ID.out.accession_id
    .map { it[1] }
    .splitCsv(header: false, sep:",")
    .map{ [it[0], it[1], it[2]] }
    .combine(CAT_ILLUMINA_FASTQ.out.reads.map { [it[0].id, it[1]] }, by: 0)
    .set { ch_sample_segment } 

  // Pull segment reference sequence for each sample
  SEQTK_SEQ(ch_sample_segment, ZSTD_DECOMPRESS_FASTA.out.file)
  ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)

  // Map reads against segment reference sequences using Minimap2
  MINIMAP2(SEQTK_SEQ.out.sample_info)
  ch_versions = ch_versions.mix(MINIMAP2.out.versions)

  // Generate coverage and stats
  MOSDEPTH_GENOME(MINIMAP2.out.alignment)
  ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions)

  FREEBAYES(MINIMAP2.out.alignment)
  ch_versions = ch_versions.mix(FREEBAYES.out.versions)

  BCF_FILTER_FREEBAYES(FREEBAYES.out.vcf, params.major_allele_fraction, params.minor_allele_fraction)
  ch_versions = ch_versions.mix(BCF_FILTER_FREEBAYES.out.versions)

  BCFTOOLS_STATS(BCF_FILTER_FREEBAYES.out.vcf)
  ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

  BCF_FILTER_FREEBAYES.out.vcf
    .combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and accession_id
    .set { ch_bcf_consensus } // ch_bcf_consensus: [sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

  COVERAGE_PLOT(ch_bcf_consensus, params.low_coverage)
  ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions)

  // Generate consensus sequences 
  BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage, params.major_allele_fraction)
  ch_versions = ch_versions.mix(BCF_CONSENSUS.out.versions)

  BCF_CONSENSUS.out.fasta
    .groupTuple(by: 0)
    .set { ch_final_consensus }

  CAT_CONSENSUS(ch_final_consensus)
  ch_versions = ch_versions.mix(CAT_CONSENSUS.out.versions)

  VADR_BCFTOOLS(CAT_CONSENSUS.out.consensus_fasta, SETUP_FLU_VADR_MODEL.out)
  ch_versions = ch_versions.mix(VADR_BCFTOOLS.out.versions)
  VADR_BCFTOOLS.out.feature_table
    .combine(VADR_BCFTOOLS.out.pass_fasta, by: 0)
    .set { ch_pre_table2asn }
  VADR_SUMMARIZE_ISSUES_BCFTOOLS(VADR_BCFTOOLS.out.vadr_outdir.map { [it[1]] }.collect())
  PRE_TABLE2ASN_BCFTOOLS(ch_pre_table2asn)
  ch_versions = ch_versions.mix(PRE_TABLE2ASN_BCFTOOLS.out.versions)
  TABLE2ASN_BCFTOOLS(PRE_TABLE2ASN_BCFTOOLS.out.table2asn_input)
  ch_versions = ch_versions.mix(TABLE2ASN_BCFTOOLS.out.versions)
  POST_TABLE2ASN_BCFTOOLS(TABLE2ASN_BCFTOOLS.out.genbank)
  ch_versions = ch_versions.mix(POST_TABLE2ASN_BCFTOOLS.out.versions)

  CAT_CONSENSUS.out.fasta
    .map { [[id:it[0]], it[1]] }
    .set { ch_cat_consensus }
	
  CAT_CONSENSUS.out.consensus_fasta
    .map { [it[0], it[1]] }
    .set { ch_cat_consensus_fasta }

  // Pass consensus sequences to GENOFLU
  GENOFLU(ch_cat_consensus_fasta)
  ch_versions = ch_versions.mix(GENOFLU.out.versions)
  
  CLEAVAGE_SITE(POST_TABLE2ASN_BCFTOOLS.out.cds_aa_fasta)
  ch_versions = ch_versions.mix(CLEAVAGE_SITE.out.versions)

  BLAST_BLASTN_CONSENSUS(ch_cat_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS.out.versions)

  ch_blastn_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
  ch_vadr_outdir_bcftools = VADR_BCFTOOLS.out.vadr_outdir.map { [it[1]] }.collect()
  SUBTYPING_REPORT_BCF_CONSENSUS(
    ZSTD_DECOMPRESS_CSV.out.file, 
    ch_blastn_consensus,
    ch_vadr_outdir_bcftools,
    CHECK_SAMPLE_SHEET.out
  )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT_BCF_CONSENSUS.out.versions)

  if (!params.skip_flumut || !params.skip_genin2) {
    PREP_FLUMUT_FASTA(CAT_CONSENSUS.out.consensus_fasta.collect({ it[1] }))
  }
  if (!params.skip_flumut) {
    FLUMUT(PREP_FLUMUT_FASTA.out.fasta)
    ch_versions = ch_versions.mix(FLUMUT.out.versions)
  }
  if (!params.skip_genin2) {
    GENIN2(PREP_FLUMUT_FASTA.out.fasta)
    ch_versions = ch_versions.mix(GENIN2.out.versions)
  }
  if (!params.skip_nextclade) {
    NEXTCLADE(
      ch_cat_consensus_fasta,
      params.nextclade_datasets_csv
    )
    ch_versions = ch_versions.mix(NEXTCLADE.out.versions)
  }

  workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
  ch_workflow_summary = Channel.value(workflow_summary)
  ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml")

  MQC_VERSIONS_TABLE(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

  // MultiQC
  ch_workflow_summary = Channel.value(Schema.params_summary_multiqc(workflow, summary_params))
  ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml")

  MULTIQC(
      ch_multiqc_config,
      MINIMAP2.out.stats.collect().ifEmpty([]),
      MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
      BCFTOOLS_STATS.out.stats.collect().ifEmpty([]),
      MQC_VERSIONS_TABLE.out.mqc_yml.collect(),
      ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
  )
}
