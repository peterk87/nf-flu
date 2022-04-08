// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LONGSHOT {
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/vcf/$sample_name",
          pattern: "*.vcf",
          mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::longshot=0.4.1 bioconda::samtools=1.15' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-e674eb3485b96f6762895eb31d74e2230db5d5bc:f69acc0b50072ac90f3d313a9d6c23487ee55ec3-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-e674eb3485b96f6762895eb31d74e2230db5d5bc:f69acc0b50072ac90f3d313a9d6c23487ee55ec3-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(bam), path(depths), path(medaka_vcf)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(depths), path(longshot_vcf), emit: vcf

    script:
    longshot_vcf = "${sample_name}.Segment_${segment}.${id}.longshot.vcf"
    """
    samtools faidx $fasta
    samtools index $bam
    longshot -P 0 -F -A --no_haps \\
        --potential_variants $medaka_vcf \\
        --bam $bam \\
        --ref $fasta \\
        --out $longshot_vcf
    """
}
