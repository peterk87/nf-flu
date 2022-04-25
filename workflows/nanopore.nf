
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_BARCODE_COUNT_FAIL } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_BARCODE_COUNT_PASS } from '../modules/local/multiqc_tsv_from_list'
include { PREPARE_NCBI_ACCESSION_ID                               } from '../modules/local/prepare_ncbi_accession_id'
include { IRMA                                                    } from '../modules/local/irma'
include { SUBTYPING_REPORT  as  SUBTYPING_REPORT_IRMA_CONSENSUS;
          SUBTYPING_REPORT  as  SUBTYPING_REPORT_BCF_CONSENSUS    } from '../modules/local/subtyping_report'
include { COVERAGE_PLOT                                           } from '../modules/local/coverage_plot'
include { VCF_FILTER_FRAMESHIFT                                   } from '../modules/local/vcf_filter_frameshift'
include { MEDAKA                                                  } from '../modules/local/medaka'
include { MINIMAP2                                                } from '../modules/local/minimap2'
include { BCF_FILTER as BCF_FILTER_CLAIR3;
          BCF_FILTER as BCF_FILTER_MEDAKA;
          BCF_CONSENSUS                                           } from '../modules/local/bcftools'
include { CLAIR3                                                  } from '../modules/local/clair3'
include { MOSDEPTH_GENOME                                         } from '../modules/local/mosdepth'
include { CAT_FASTQ;
          GUNZIP as GUNZIP_FLU_FASTA;
          CAT_DB; CAT_CONSENSUS                                   } from '../modules/local/misc'
include { BLAST_BLASTDBCMD                                        } from '../modules/local/pull_references'
include { CHECK_SAMPLE_SHEET                                      } from '../modules/local/check_sample_sheet'
include { REF_FASTA_CHECK                                         } from '../modules/local/ref_fasta_check'

include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI_NO_PARSEID  } from '../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_REFDB            } from '../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI_PARSEID     } from '../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_BLASTN_NCBI                       } from '../modules/nf-core/modules/blast/blastn/main'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS                  } from '../modules/nf-core/modules/blast/blastn/main'

if (params.input) { ch_input= file(params.input) }

def pass_barcode_reads = [:]
def fail_barcode_reads = [:]
ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)
def irma_module = 'FLU-minion'
if (params.irma_module) {
    irma_module = params.irma_module
}

workflow NANOPORE {
    ch_versions = Channel.empty()
    ch_input_ref = Channel.empty()

    Channel.fromPath(params.input, checkIfExists: true)
    | CHECK_SAMPLE_SHEET
    | splitCsv(header: ['sample', 'barcode'], sep: ',', skip: 1) \
    | map { it ->
        def count = 0
        def fastqFiles = file(it.barcode).listFiles()
        for (x in fastqFiles) {
            if (x.isFile() && x.toString().contains('.fastq')) {
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
    MULTIQC_TSV_BARCODE_COUNT_FAIL (
        ch_pass_fail_barcode_count.fail.collect(),
        ['Sample', 'Barcode count'],
        'fail_barcode_count_samples'
    )
    // Report samples which have reads count  > min_barcode_reads
    MULTIQC_TSV_BARCODE_COUNT_PASS (
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

    ch_input_ref_db = GUNZIP_FLU_FASTA.out.gunzip

    if (params.ref_db){
        ref_fasta_file = file(params.ref_db, type: 'file')
        REF_FASTA_CHECK(ref_fasta_file)
        CAT_DB(GUNZIP_FLU_FASTA.out.gunzip, REF_FASTA_CHECK.out.fasta)
        ch_input_ref_db = CAT_DB.out.fasta
        //BLAST_MAKEBLASTDB_REFDB(REF_FASTA_CHECK.out.fasta)
    }

    BLAST_MAKEBLASTDB_NCBI_NO_PARSEID(ch_input_ref_db)

    // Make another blast db with parse_id option for pulling reference based on Accession ID
    BLAST_MAKEBLASTDB_NCBI_PARSEID(ch_input_ref_db)

    CAT_FASTQ(ch_fastq_dirs)

    CAT_FASTQ.out.reads
    | map { it ->
        return [it[0].id, it[1]]
    }
    | set {ch_aggregate_reads}

    if (params.ref_sequences){
        //ref_sequences in format: ref_id, segment_no, ref_fasta
        Channel.fromPath(params.ref_sequences, checkIfExists: true) \
        | splitCsv(header: false) \
        | map { it ->
            // 0: ref ID, 1: segment, 2: path ref_fasta
            return [it[0], it[1], it[2]]
        } \
        | combine (ch_aggregate_reads) \
        | set {ch_input_ref_combine }

        ch_input_ref_combine
        | map { it ->
            // get [sample, segment, ref_id, ref_fasta, reads]
            return [it[3], it[1], it[0], it[2], it[4]]
        } | set { ch_input_ref }
    }

    // IRMA to generate amended consensus sequences
    IRMA(CAT_FASTQ.out.reads, irma_module)
    //ch_versions.mix(IRMA.out.version)
    // Find the top map sequences against ncbi database
    BLAST_BLASTN_NCBI(IRMA.out.consensus, BLAST_MAKEBLASTDB_NCBI_NO_PARSEID.out.db)
    //ch_versions.mix(BLAST_BLASTN.out.versions)

    //Generate suptype prediction report
    /* Does not need for now
    ch_blast = BLAST_BLASTN_NCBI.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_IRMA_CONSENSUS(ch_influenza_metadata, ch_blast)
    */

    // Prepare top ncbi accession id for each segment of each sample sample (id which has top bitscore)
    PREPARE_NCBI_ACCESSION_ID(BLAST_BLASTN_NCBI.out.txt, ch_influenza_metadata)

    PREPARE_NCBI_ACCESSION_ID.out.accession_id
    | map {it[1]} | splitCsv(header: false, sep:",")
    | map{ it ->
        // 0: sample_name, 1: segment, 2: ref_ncbi_accession_id
        return [it[0],it[1],it[2]]
    }
    | combine(ch_aggregate_reads, by: 0)
    | set {ch_sample_segment} // ch_sample_segment: [sample_name, segment, id, reads]

    //Pull segment reference sequence for each sample
    BLAST_BLASTDBCMD(ch_sample_segment, BLAST_MAKEBLASTDB_NCBI_PARSEID.out.db)

    // Map reads against segment reference sequences
    ch_mapping = BLAST_BLASTDBCMD.out.fasta.concat(ch_input_ref)
    MINIMAP2(ch_mapping)
    //ch_versions.mix(MINIMAP2.out.version)

    MOSDEPTH_GENOME(MINIMAP2.out.alignment)
    //ch_versions.mix(MOSDEPTH_GENOME.out.version)

    // Variants calling
    MINIMAP2.out.alignment
    | map { it->
        // [0: sample_name, segment, 2: id, 3: fasta, 4: path('*.{bam,bam.bai}'), 5:depths]
        return [it[0], it[1], it[2], it[3], it[4], it[5]]
    } | set {ch_variant_calling}

    if (!params.skip_clair3){
        CLAIR3(ch_variant_calling)
        BCF_FILTER_CLAIR3(CLAIR3.out.vcf, params.major_allele_fraction)
        ch_vcf_filter = BCF_FILTER_CLAIR3.out.vcf
    } else {
        MEDAKA(ch_variant_calling)
        BCF_FILTER_MEDAKA(MEDAKA.out.vcf, params.major_allele_fraction)
        ch_vcf_filter = BCF_FILTER_MEDAKA.out.vcf
    }

    VCF_FILTER_FRAMESHIFT(ch_vcf_filter)

    COVERAGE_PLOT (VCF_FILTER_FRAMESHIFT.out.vcf, params.low_coverage)

    VCF_FILTER_FRAMESHIFT.out.vcf
    | map { it ->
        //[0: sample_name, 1: segment, 2: id, 3: fasta, 4: depths, 5: filt_vcf]
        return [it[0], it[1], it[2], it[3], it[5]] // no need 4: depths
    }
    | combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and accession_id
    | set { ch_bcf_consensus } //[sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

    // Generate consensus sequences
    BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)
    BCF_CONSENSUS.out.fasta
    | groupTuple(by: 0)
    | set { ch_final_consensus }

    CAT_CONSENSUS(ch_final_consensus)
    CAT_CONSENSUS.out.fasta
    | map {it ->
        return [[id:it[0]], it[1]]
    }
    | set { ch_blatn_consensus }
    BLAST_BLASTN_CONSENSUS(ch_blatn_consensus, BLAST_MAKEBLASTDB_NCBI_NO_PARSEID.out.db)
    ch_blast_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_BCF_CONSENSUS(ch_influenza_metadata, ch_blast_consensus)
}