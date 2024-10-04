#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV        } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV        } from '../modules/local/multiqc_tsv_from_list'
include { PULL_TOP_REF_ID                                     } from '../modules/local/pull_top_ref_id'
include { IRMA                                                } from '../modules/local/irma'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_IRMA_CONSENSUS } from '../modules/local/subtyping_report'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_BCF_CONSENSUS  } from '../modules/local/subtyping_report'
include { COVERAGE_PLOT                                       } from '../modules/local/coverage_plot'
include { BLASTN_REPORT                                       } from '../modules/local/blastn_report'
include { VCF_FILTER_FRAMESHIFT                               } from '../modules/local/vcf_filter_frameshift'
include { MEDAKA                                              } from '../modules/local/medaka'
include { MINIMAP2                                            } from '../modules/local/minimap2'
include { BCF_FILTER as BCF_FILTER_CLAIR3                     } from '../modules/local/bcftools'
include { BCF_FILTER as BCF_FILTER_MEDAKA                     } from '../modules/local/bcftools'
include { BCF_CONSENSUS; BCFTOOLS_STATS                       } from '../modules/local/bcftools'
include { CLAIR3                                              } from '../modules/local/clair3'
include { MOSDEPTH_GENOME                                     } from '../modules/local/mosdepth'
include { CAT_NANOPORE_FASTQ                                  } from '../modules/local/misc'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_FASTA; ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_CSV } from '../modules/local/zstd_decompress'
include { CAT_DB                                              } from '../modules/local/misc'
include { CAT_CONSENSUS                                       } from '../modules/local/misc'
include { SEQTK_SEQ                                           } from '../modules/local/seqtk_seq'
include { CHECK_SAMPLE_SHEET                                  } from '../modules/local/check_sample_sheet'
include { CHECK_REF_FASTA                                     } from '../modules/local/check_ref_fasta'
// using modified BLAST_MAKEBLASTDB from nf-core/modules to only move/publish BLAST DB files
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI         } from '../modules/local/blast_makeblastdb'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_REFDB        } from '../modules/local/blast_makeblastdb'
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                   } from '../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS              } from '../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS_REF_DB       } from '../modules/local/blastn'
include { SETUP_FLU_VADR_MODEL; VADR; VADR_SUMMARIZE_ISSUES   } from '../modules/local/vadr'
include { PRE_TABLE2ASN; TABLE2ASN; POST_TABLE2ASN            } from '../modules/local/table2asn'
include { MULTIQC                                             } from '../modules/local/multiqc'
include { MQC_VERSIONS_TABLE } from '../modules/local/mqc_versions_table'

def pass_sample_reads = [:]
def fail_sample_reads = [:]
ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)
ch_vadr_model_targz = file(params.vadr_model_targz)
if (params.clair3_user_variant_model) {
  ch_user_clair3_model = file(params.clair3_user_variant_model, checkIfExists: true)
}
def irma_module = 'FLU-minion'
if (params.irma_module) {
    irma_module = params.irma_module
}
def json_schema = "$projectDir/nextflow_schema.json"
def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)

workflow NANOPORE {
  ch_versions = Channel.empty()

  ch_input = CHECK_SAMPLE_SHEET(Channel.fromPath(params.input, checkIfExists: true))

  ch_input.splitCsv(header: ['sample', 'reads'], sep: ',', skip: 1)
    // "reads" can be path to file or directory
    .map { [it.sample, it.reads] }
    // group by sample name to later merge all reads for that sample
    .groupTuple(by: 0)
    // collect all uncompressed and compressed FASTQ reads into 2 lists
    // and count number of reads for sample
    .map { sample, reads ->
      // uncompressed FASTQ list
      def fq = []
      // compressed FASTQ list
      def fqgz = []
      // read count
      def count = 0
      for (f in reads) {
        f = file(f)
        if (f.isFile() && f.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
          if (f.getName() ==~ /.*\.gz/) {
            fqgz << f
          } else {
            fq << f
          }
          continue
        }
        if (f.isDirectory()) {
          // only look for FQ reads in first level of directory
          for (x in f.listFiles()) {
            if (x.isFile() && x.getName() ==~ /.*\.(fastq|fq)(\.gz)?/) {
              if (x.getName() ==~ /.*\.gz/) {
                fqgz << x
              } else {
                fq << x
              }
            }
          }
        }
      }
      for (x in fq) {
        count += x.countFastq()
      }
      for (x in fqgz) {
        count += x.countFastq()
      }
      return [ sample, fqgz, fq, count ]
    }
    .set { ch_input_sorted }

  ch_input_sorted
    .branch { sample, fqgz, fq, count  ->
      pass: count >= params.min_sample_reads
        pass_sample_reads[sample] = count
        return [ "$sample\t$count" ]
      fail: count < params.min_sample_reads
        fail_sample_reads[sample] = count
        return [ "$sample\t$count" ]
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

  // Keep samples which have reads count  > min_sample_reads for downstream analysis
  // Re-arrange channels to have meta map of information for sample
  ch_input_sorted
    .filter { it[-1] >= params.min_sample_reads }
    .map { sample, fqgz, fq, count -> [ [id: sample], fqgz, fq ] }
    .set { ch_reads }

  ZSTD_DECOMPRESS_FASTA(ch_influenza_db_fasta, "influenza.fasta")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_FASTA.out.versions)
  ZSTD_DECOMPRESS_CSV(ch_influenza_metadata, "influenza.csv")
  ch_versions = ch_versions.mix(ZSTD_DECOMPRESS_CSV.out.versions)

  ch_input_ref_db = ZSTD_DECOMPRESS_FASTA.out.file

  if (params.ref_db){
    ch_ref_fasta = file(params.ref_db, type: 'file')
    CHECK_REF_FASTA(ch_ref_fasta)
    ch_versions = ch_versions.mix(CHECK_REF_FASTA.out.versions)
    CAT_DB(ZSTD_DECOMPRESS_FASTA.out.file, CHECK_REF_FASTA.out.fasta)
    ch_input_ref_db = CAT_DB.out.fasta
  }

  BLAST_MAKEBLASTDB_NCBI(ch_input_ref_db)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_NCBI.out.versions)

  SETUP_FLU_VADR_MODEL(ch_vadr_model_targz)

  CAT_NANOPORE_FASTQ(ch_reads)

  // IRMA to generate amended consensus sequences
  IRMA(CAT_NANOPORE_FASTQ.out.reads, irma_module)
  ch_versions = ch_versions.mix(IRMA.out.versions)
  // Find the top map sequences against ncbi database
  BLAST_BLASTN_IRMA(IRMA.out.majority_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN_IRMA.out.versions)

  // Generate suptype prediction report from IRMA results
  if (!params.skip_irma_subtyping_report){
    ch_blast_irma = BLAST_BLASTN_IRMA.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_IRMA_CONSENSUS(
      ZSTD_DECOMPRESS_CSV.out.file,
      ch_blast_irma,
      CHECK_SAMPLE_SHEET.out
    )
  }

  // Prepare top ncbi accession id for each segment of each sample sample (id which has top bitscore)
  PULL_TOP_REF_ID(BLAST_BLASTN_IRMA.out.txt, ZSTD_DECOMPRESS_CSV.out.file)
  ch_versions = ch_versions.mix(PULL_TOP_REF_ID.out.versions)

  PULL_TOP_REF_ID.out.accession_id
    .map { it[1] }
    .splitCsv(header: false, sep:",")
    // 0: sample_name, 1: segment, 2: ref_ncbi_accession_id, 3: ref_sequence_name
    .map{ [it[0], it[1], it[2]] }
    .combine(CAT_NANOPORE_FASTQ.out.reads.map { [it[0].id, it[1]] }, by: 0)
    .set { ch_sample_segment } // ch_sample_segment: [sample_name, segment, id, reads]

  // Pull segment reference sequence for each sample
  SEQTK_SEQ(ch_sample_segment, ch_input_ref_db)
  ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)

  // Map reads against segment reference sequences
  MINIMAP2(SEQTK_SEQ.out.sample_info)
  ch_versions = ch_versions.mix(MINIMAP2.out.versions)

  MOSDEPTH_GENOME(MINIMAP2.out.alignment)
  ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions)

  // Variants calling
  if (params.variant_caller == 'clair3'){
    if (params.clair3_user_variant_model) {
      CLAIR3(
        MINIMAP2.out.alignment,
        ch_user_clair3_model
      )
    } else {
      CLAIR3(MINIMAP2.out.alignment, [])
    }

    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    BCF_FILTER_CLAIR3(CLAIR3.out.vcf, params.major_allele_fraction, params.minor_allele_fraction)
    ch_versions = ch_versions.mix(BCF_FILTER_CLAIR3.out.versions)
    ch_vcf_filter = BCF_FILTER_CLAIR3.out.vcf
  } else if (params.variant_caller == 'medaka') {
    MEDAKA(MINIMAP2.out.alignment)
    ch_versions = ch_versions.mix(MEDAKA.out.versions)

    BCF_FILTER_MEDAKA(MEDAKA.out.vcf, params.major_allele_fraction, params.minor_allele_fraction)
    ch_versions = ch_versions.mix(BCF_FILTER_MEDAKA.out.versions)
    ch_vcf_filter = BCF_FILTER_MEDAKA.out.vcf
  }

  VCF_FILTER_FRAMESHIFT(ch_vcf_filter)
  ch_versions = ch_versions.mix(VCF_FILTER_FRAMESHIFT.out.versions)

  BCFTOOLS_STATS(VCF_FILTER_FRAMESHIFT.out.vcf)
  ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

  VCF_FILTER_FRAMESHIFT.out.vcf
    .combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and accession_id
    .set { ch_bcf_consensus } // ch_bcf_consensus: [sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

  COVERAGE_PLOT(ch_bcf_consensus, params.low_coverage)
  ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions)

  // Generate consensus sequences
  BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)
  ch_versions = ch_versions.mix(BCF_CONSENSUS.out.versions)

  BCF_CONSENSUS.out.fasta
    .groupTuple(by: 0)
    .set { ch_final_consensus }

  CAT_CONSENSUS(ch_final_consensus)
  ch_versions = ch_versions.mix(CAT_CONSENSUS.out.versions)

  VADR(CAT_CONSENSUS.out.consensus_fasta, SETUP_FLU_VADR_MODEL.out)
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

  CAT_CONSENSUS.out.fasta
    .map { [[id:it[0]], it[1]] }
    .set { ch_cat_consensus }

  BLAST_BLASTN_CONSENSUS(ch_cat_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS.out.versions)

  ch_blastn_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
  SUBTYPING_REPORT_BCF_CONSENSUS(
    ZSTD_DECOMPRESS_CSV.out.file, 
    ch_blastn_consensus,
    CHECK_SAMPLE_SHEET.out
  )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT_BCF_CONSENSUS.out.versions)

  if (params.ref_db){
    BLAST_MAKEBLASTDB_REFDB(CHECK_REF_FASTA.out.fasta)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_REFDB.out.versions)

    BLAST_BLASTN_CONSENSUS_REF_DB(ch_cat_consensus, BLAST_MAKEBLASTDB_REFDB.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS_REF_DB.out.versions)

    BLASTN_REPORT(BLAST_BLASTN_CONSENSUS_REF_DB.out.txt)
    ch_versions = ch_versions.mix(BLASTN_REPORT.out.versions)
  }

  workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
  ch_workflow_summary = Channel.value(workflow_summary)
  ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml")

  MQC_VERSIONS_TABLE(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

  MULTIQC(
      ch_multiqc_config,
      MINIMAP2.out.stats.collect().ifEmpty([]),
      MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
      BCFTOOLS_STATS.out.stats.collect().ifEmpty([]),
      MQC_VERSIONS_TABLE.out.mqc_yml.collect(),
      ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
  )
}
