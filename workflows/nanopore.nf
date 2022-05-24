
include { MULTIQC_TSV_FROM_LIST as BARCODE_COUNT_FAIL             } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as BARCODE_COUNT_PASS             } from '../modules/local/multiqc_tsv_from_list'
include { PULL_TOP_REF_ID                                         } from '../modules/local/pull_top_ref_id'
include { IRMA                                                    } from '../modules/local/irma'
include { SUBTYPING_REPORT  as  SUBTYPING_REPORT_IRMA_CONSENSUS;
          SUBTYPING_REPORT  as  SUBTYPING_REPORT_BCF_CONSENSUS    } from '../modules/local/subtyping_report'
include { COVERAGE_PLOT                                           } from '../modules/local/coverage_plot'
include { BLASTN_REPORT                                           } from '../modules/local/blastn_report'
include { VCF_FILTER_FRAMESHIFT                                   } from '../modules/local/vcf_filter_frameshift'
include { MEDAKA                                                  } from '../modules/local/medaka'
include { MINIMAP2                                                } from '../modules/local/minimap2'
include { BCF_FILTER as BCF_FILTER_CLAIR3;
          BCF_FILTER as BCF_FILTER_MEDAKA;
          BCF_CONSENSUS; BCFTOOLS_STATS                           } from '../modules/local/bcftools'
include { CLAIR3                                                  } from '../modules/local/clair3'
include { MOSDEPTH_GENOME                                         } from '../modules/local/mosdepth'
include { CAT_FASTQ;
          GUNZIP as GUNZIP_FLU_FASTA;
          CAT_DB; CAT_CONSENSUS                                   } from '../modules/local/misc'
include { SEQTK_SEQ                                               } from '../modules/local/seqtk_seq'
include { CHECK_SAMPLE_SHEET                                      } from '../modules/local/check_sample_sheet'
include { CHECK_REF_FASTA                                         } from '../modules/local/check_ref_fasta'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI             } from '../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_REFDB            } from '../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                       } from '../modules/nf-core/modules/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS                  } from '../modules/nf-core/modules/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS_REF_DB           } from '../modules/nf-core/modules/blast/blastn/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS  as SOFTWARE_VERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC as MULTIQC_NANOPORE                             } from '../modules/local/multiqc'

def pass_barcode_reads = [:]
def fail_barcode_reads = [:]
ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)
def irma_module = 'FLU-minion'
if (params.irma_module) {
    irma_module = params.irma_module
}
def json_schema = "$projectDir/nextflow_schema.json"
def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)

workflow NANOPORE {
    ch_versions = Channel.empty()

    Channel.fromPath(params.input, checkIfExists: true)
    | CHECK_SAMPLE_SHEET
    | splitCsv(header: ['sample', 'barcode'], sep: ',', skip: 1) \
    | map { it ->
        def count = 0
        def fastqFiles = file(it.barcode).listFiles()
        for (x in fastqFiles) {
            if (x.isFile() && (x.toString().contains('.fastq') || x.toString().contains('.fq'))) {
                count += x.countFastq()
            }
        }
        return [it.sample , it.barcode, count]
    } | set { ch_fastq_dirs }

    ch_fastq_dirs
    | branch { sample, dir, count  ->
        pass: count > params.min_barcode_reads
            pass_barcode_reads[sample] = count
            return [ "$sample\t$count" ]
        fail: count <= params.min_barcode_reads
            fail_barcode_reads[sample] = count
            return [ "$sample\t$count" ]
    }
    | set { ch_pass_fail_barcode_count }

    // Report samples which have reads count  <= min_barcode_reads
    BARCODE_COUNT_FAIL (
        ch_pass_fail_barcode_count.fail.collect(),
        ['Sample', 'Barcode count'],
        'fail_barcode_count_samples'
    )
    // Report samples which have reads count  > min_barcode_reads
    BARCODE_COUNT_PASS (
        ch_pass_fail_barcode_count.pass.collect(),
        ['Sample', 'Barcode count'],
        'pass_barcode_count_samples'
    )

    // Keep samples which have reads count  > min_barcode_reads for downstream analysis
    // Re-arrange channels to have meta map of information for sample
    ch_fastq_dirs
    | filter { it[-1] > params.min_barcode_reads }
    | map { sample, dir, count -> [ [ id: sample], dir ] }
    | set { ch_fastq_dirs }

    GUNZIP_FLU_FASTA(ch_influenza_db_fasta)
    ch_versions = ch_versions.mix(GUNZIP_FLU_FASTA.out.versions)

    ch_input_ref_db = GUNZIP_FLU_FASTA.out.gunzip

    if (params.ref_db){
        ref_fasta_file = file(params.ref_db, type: 'file')
        CHECK_REF_FASTA(ref_fasta_file)
        ch_versions = ch_versions.mix(CHECK_REF_FASTA.out.versions)
        CAT_DB(GUNZIP_FLU_FASTA.out.gunzip, CHECK_REF_FASTA.out.fasta)
        ch_input_ref_db = CAT_DB.out.fasta
    }

    BLAST_MAKEBLASTDB_NCBI(ch_input_ref_db)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_NCBI.out.versions)

    CAT_FASTQ(ch_fastq_dirs)
    CAT_FASTQ.out.reads
    | map { it ->
        return [it[0].id, it[1]]
    }
    | set {ch_aggregate_reads}

    // IRMA to generate amended consensus sequences
    IRMA(CAT_FASTQ.out.reads, irma_module)
    ch_versions = ch_versions.mix(IRMA.out.versions)
    // Find the top map sequences against ncbi database
    BLAST_BLASTN_IRMA(IRMA.out.consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN_IRMA.out.versions)

    //Generate suptype prediction report
    if (!params.skip_irma_subtyping_report){
        ch_blast_irma = BLAST_BLASTN_IRMA.out.txt.collect({ it[1] })
        SUBTYPING_REPORT_IRMA_CONSENSUS(ch_influenza_metadata, ch_blast_irma)
    }

    // Prepare top ncbi accession id for each segment of each sample sample (id which has top bitscore)
    PULL_TOP_REF_ID(BLAST_BLASTN_IRMA.out.txt, ch_influenza_metadata)
    ch_versions = ch_versions.mix(PULL_TOP_REF_ID.out.versions)

    PULL_TOP_REF_ID.out.accession_id
    | map {it[1]} | splitCsv(header: false, sep:",")
    | map{ it ->
        // 0: sample_name, 1: segment, 2: ref_ncbi_accession_id, 3: ref_sequence_name
        return [it[0],it[1],it[2]]
    }
    | combine(ch_aggregate_reads, by: 0)
    | set {ch_sample_segment} // ch_sample_segment: [sample_name, segment, id, reads]

    //Pull segment reference sequence for each sample
    SEQTK_SEQ(ch_sample_segment, ch_input_ref_db)
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)

    // Map reads against segment reference sequences
    ch_mapping = SEQTK_SEQ.out.sample_info

    MINIMAP2(ch_mapping)
    ch_versions = ch_versions.mix(MINIMAP2.out.versions)

    MOSDEPTH_GENOME(MINIMAP2.out.alignment)
    ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions)

    // Variants calling
    MINIMAP2.out.alignment
    | map { it->
        // [0: sample_name, segment, 2: id, 3: fasta, 4: path('*.{bam,bam.bai}'), 5:depths]
        return [it[0], it[1], it[2], it[3], it[4], it[5]]
    } | set {ch_variant_calling}

    if (params.variant_caller == 'clair3'){
        CLAIR3(ch_variant_calling)
        ch_versions = ch_versions.mix(CLAIR3.out.versions)

        BCF_FILTER_CLAIR3(CLAIR3.out.vcf, params.major_allele_fraction)
        ch_versions = ch_versions.mix(BCF_FILTER_CLAIR3.out.versions)
        ch_vcf_filter = BCF_FILTER_CLAIR3.out.vcf
    } else if (params.variant_caller == 'medaka') {
        MEDAKA(ch_variant_calling)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)

        BCF_FILTER_MEDAKA(MEDAKA.out.vcf, params.major_allele_fraction)
        ch_versions = ch_versions.mix(BCF_FILTER_MEDAKA.out.versions)
        ch_vcf_filter = BCF_FILTER_MEDAKA.out.vcf
    }

    VCF_FILTER_FRAMESHIFT(ch_vcf_filter)
    ch_versions = ch_versions.mix(VCF_FILTER_FRAMESHIFT.out.versions)

    BCFTOOLS_STATS(VCF_FILTER_FRAMESHIFT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    COVERAGE_PLOT(VCF_FILTER_FRAMESHIFT.out.vcf, params.low_coverage)
    ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions)

    VCF_FILTER_FRAMESHIFT.out.vcf
    | map { it ->
        //[0: sample_name, 1: segment, 2: id, 3: fasta, 4: depths, 5: filt_vcf]
        return [it[0], it[1], it[2], it[3], it[5]] // no need 4: depths
    }
    | combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and accession_id
    | set { ch_bcf_consensus } //[sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

    // Generate consensus sequences
    BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)
    ch_versions = ch_versions.mix(BCF_CONSENSUS.out.versions)

    BCF_CONSENSUS.out.fasta
    | groupTuple(by: 0)
    | set { ch_final_consensus }

    CAT_CONSENSUS(ch_final_consensus)
    ch_versions = ch_versions.mix(CAT_CONSENSUS.out.versions)

    CAT_CONSENSUS.out.fasta
    | map {it ->
        return [[id:it[0]], it[1]]
    }
    | set { ch_cat_consensus }

    BLAST_BLASTN_CONSENSUS(ch_cat_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS.out.versions)

    ch_blastn_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_BCF_CONSENSUS(ch_influenza_metadata, ch_blastn_consensus)
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
    SOFTWARE_VERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    MULTIQC_NANOPORE(
        ch_multiqc_config,
        MINIMAP2.out.stats.collect().ifEmpty([]),
        MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
        BCFTOOLS_STATS.out.stats.collect().ifEmpty([]),
        SOFTWARE_VERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    )
}