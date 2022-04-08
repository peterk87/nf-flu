// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MEDAKA{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/vcf/$sample_name",
        pattern: "*.vcf",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::medaka=1.6.0 bioconda::samtools=1.15' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-1d9673528d7dd295349af77d8ecfaad41467ae6d:5ff5398dd011d6ff444bea563bd6dd4fa8cc33b2-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-1d9673528d7dd295349af77d8ecfaad41467ae6d:5ff5398dd011d6ff444bea563bd6dd4fa8cc33b2-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(bam), path(depths), path(flagstat), path(idxstats)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(bam), path(depths), path(vcf), emit: stat

    script:
    vcf = "${sample_name}.Segment_${segment}.${id}.medaka.vcf"
    """
    samtools index $bam
    medaka consensus --chunk_len 800 --chunk_ovlp 400 $bam ${bam}.hdf
    medaka variant $fasta ${bam}.hdf $vcf
    """
}
