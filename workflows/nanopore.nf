include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_NO_SAMPLE_NAME     } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_NO_BARCODES        } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_BARCODE_COUNT      } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_GUPPYPLEX_COUNT    } from '../modules/local/multiqc_tsv_from_list'
include { PREPARE_NCBI_ACCESSION_ID                               } from '../modules/local/prepare_ncbi_accession_id'
include { IRMA                                                    } from '../modules/local/irma'
include { SUBTYPING_REPORT                                        } from '../modules/local/subtyping_report'
include { COVERAGE_PLOT                                           } from '../modules/local/coverage_plot'
include { VCF_FILTER_FRAMESHIFT                                   } from '../modules/local/vcf_filter_frameshift'
include { MEDAKA                                                  } from '../modules/local/medaka'
include { MINIMAP2                                                } from '../modules/local/minimap2'
include { BCF_FILTER; BCF_CONSENSUS                               } from '../modules/local/bcftools'
include { MOSDEPTH_GENOME                                         } from '../modules/local/mosdepth'
include { CAT_FASTQ; GUNZIP as GUNZIP_FLU_FASTA                   } from '../modules/local/misc'
include { BLAST_BLASTDBCMD                                        } from '../modules/local/pull_references'
include { INPUT_CHECK                                             } from '../subworkflows/local/input_check'

include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NO_PARSEID       } from '../modules/nf-core/modules/blast/makeblastdb/main' //addParams( options: modules['blastn_makeblastdb'] )
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_PARSEID          } from '../modules/nf-core/modules/blast/makeblastdb/main' //addParams( options: modules['blastn_makeblastdb_parseid'] )
include { BLAST_BLASTN                                            } from '../modules/nf-core/modules/blast/blastn/main'      //addParams( options: modules['blast_blastn'] )

if (params.input) { ch_input= file(params.input) }

def pass_barcode_reads = [:]
def fail_barcode_reads = [:]
ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)

workflow NANOPORE {
    ch_versions = Channel.empty()

    barcode_dirs       = file("${params.fastq_dir}/barcode*", type: 'dir' , maxdepth: 1)
    single_barcode_dir = file("${params.fastq_dir}/*.fastq" , type: 'file', maxdepth: 1)
    ch_custom_no_sample_name_multiqc = Channel.empty()
    ch_custom_no_barcodes_multiqc    = Channel.empty()

    if (barcode_dirs) {
        Channel
            .fromPath( barcode_dirs )
            .filter( ~/.*barcode[0-9]{1,4}$/ )
            .map { dir ->
                def count = 0
                for (x in dir.listFiles()) {
                    if (x.isFile() && x.toString().contains('.fastq')) {
                        count += x.countFastq()
                    }
                }
                return [dir.baseName , dir, count]
            }
            .set { ch_fastq_dirs }

        //
        // SUBWORKFLOW: Read in samplesheet containing sample to barcode mappings
        //
        if (params.input) {
            INPUT_CHECK (
                ch_input,
                params.platform
            )
            .sample_info
            .join(ch_fastq_dirs, remainder: true)
            .set { ch_fastq_dirs }

            //
            // MODULE: Create custom content file for MultiQC to report barcodes were allocated reads >= params.min_barcode_reads but no sample name in samplesheet
            //
            ch_fastq_dirs
                .filter { it[1] == null }
                .filter { it[-1] >= params.min_barcode_reads }
                .map { it -> [ "${it[0]}\t${it[-1]}" ] }
                .set { ch_barcodes_no_sample }

            MULTIQC_TSV_NO_SAMPLE_NAME (
                ch_barcodes_no_sample.collect(),
                ['Barcode', 'Read count'],
                'fail_barcodes_no_sample'
            )
            .set { ch_custom_no_sample_name_multiqc }

            //
            // MODULE: Create custom content file for MultiQC to report samples that were in samplesheet but have no barcodes
            //
            ch_fastq_dirs
                .filter { it[-1] == null }
                .map { it -> [ "${it[1]}\t${it[0]}" ] }
                .set { ch_samples_no_barcode }

            MULTIQC_TSV_NO_BARCODES (
                ch_samples_no_barcode.collect(),
                ['Sample', 'Missing barcode'],
                'fail_no_barcode_samples'
            )
            .set { ch_custom_no_barcodes_multiqc }

            ch_fastq_dirs
                .filter { (it[1] != null)  }
                .filter { (it[-1] != null) }
                .set { ch_fastq_dirs }

        } else {
            ch_fastq_dirs
                .map { barcode, dir, count -> [ barcode, barcode, dir, count ] }
                .set { ch_fastq_dirs }
        }
    } else if (single_barcode_dir) {
        Channel
            .fromPath("${params.fastq_dir}", type: 'dir', maxDepth: 1)
            .map { it -> [ 'SAMPLE_1', 'single_barcode', it, 10000000 ] }
            .set{ ch_fastq_dirs }
    } else {
        log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files generated by guppy_barcoder or guppy_basecaller e.g. '--fastq_dir ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/fastq_pass/"
        System.exit(1)
    }

    //
    // MODULE: Create custom content file for MultiQC to report samples with reads < params.min_barcode_reads
    //
    ch_fastq_dirs
        .branch { barcode, sample, dir, count  ->
            pass: count > params.min_barcode_reads
                pass_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
            fail: count < params.min_barcode_reads
                fail_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
        }
        .set { ch_pass_fail_barcode_count }

    MULTIQC_TSV_BARCODE_COUNT (
        ch_pass_fail_barcode_count.fail.collect(),
        ['Sample', 'Barcode count'],
        'fail_barcode_count_samples'
    )

    // Re-arrange channels to have meta map of information for sample
    ch_fastq_dirs
        .filter { it[-1] > params.min_barcode_reads }
        .map { barcode, sample, dir, count -> [ [ id: sample, barcode:barcode ], dir ] }
        .set { ch_fastq_dirs }


    GUNZIP_FLU_FASTA(ch_influenza_db_fasta)

    BLAST_MAKEBLASTDB_NO_PARSEID(GUNZIP_FLU_FASTA.out.gunzip)

    // Make another blast db with parse_id option for pulling reference based on Accession ID
    BLAST_MAKEBLASTDB_PARSEID(GUNZIP_FLU_FASTA.out.gunzip)

    CAT_FASTQ(ch_fastq_dirs).set {ch_aggregate_reads}

    ch_aggregate_reads.map { it ->
        return [it[0].id, it[1]]
    }
    .set {ch_aggregate_reads}

    // IRMA to generate amended consensus sequences
    IRMA(CAT_FASTQ.out.reads)
    ch_versions.mix(IRMA.out.version)

    // Find the top map sequences against ncbi database
    BLAST_BLASTN(IRMA.out.consensus, BLAST_MAKEBLASTDB_NO_PARSEID.out.db)
    ch_versions.mix(BLAST_BLASTN.out.versions)

    //Generate suptype prediction report
    ch_blast = BLAST_BLASTN.out.txt.collect({ it[1] })
    SUBTYPING_REPORT(ch_influenza_metadata, ch_blast)

    // Prepare top ncbi accession id for each segment of each sample sample (id which has top bitscore)
    PREPARE_NCBI_ACCESSION_ID(BLAST_BLASTN.out.txt, ch_influenza_metadata)
    .set {ch_ref_accession_id}


    ch_ref_accession_id.map {it[1]}.splitCsv(sep:",")
    .map{ it ->
        return [it[0],it[1],it[2]]
    }
    .combine(ch_aggregate_reads, by: 0)
    .set {ch_sample_segment} // ch_sample_segment: [sample_name, segment, id, reads]

    //Pull segment reference sequence for each sample
    BLAST_BLASTDBCMD(ch_sample_segment, BLAST_MAKEBLASTDB_PARSEID.out.db)

    // Map reads againts segment reference sequences
    MINIMAP2(BLAST_BLASTDBCMD.out.fasta)
    ch_versions.mix(MINIMAP2.out.version)

    MOSDEPTH_GENOME(MINIMAP2.out.alignment)
    ch_versions.mix(MOSDEPTH_GENOME.out.version)

    // Variants calling
    MINIMAP2.out.alignment
    .map { it->
        // [0: sample_name, segment, 2: id, 3: fasta, 4: path('*.{bam,bam.bai}'), 5:depths]
        return [it[0], it[1], it[2], it[3], it[4], it[5]]
    }. set {ch_medaka}
    MEDAKA(ch_medaka)
    ch_versions.mix(MEDAKA.out.version)

    BCF_FILTER(MEDAKA.out.vcf, params.major_allele_fraction)
    ch_versions = ch_versions.mix(BCF_FILTER.out.version)

    VCF_FILTER_FRAMESHIFT(BCF_FILTER.out.vcf)

    COVERAGE_PLOT (VCF_FILTER_FRAMESHIFT.out.vcf, params.low_coverage)

    VCF_FILTER_FRAMESHIFT.out.vcf
    .map { it ->
        //[0: sample_name, 1: segment, 2: id, 3: fasta, 4: depths, 5: filt_vcf]
        return [it[0], it[1], it[2], it[3], it[5]] // no need 4: depths
    }
    .combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and acession_id
    .set { ch_bcf_consensus } //[sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

    // Generate consensus sequences
    BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)

}