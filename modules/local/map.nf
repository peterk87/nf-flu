// Import generic module functions
include { getSoftwareName } from './functions'

process MAP{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_high'
    publishDir "${params.outdir}/mapping/$sample_name",
         mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::minimap2=2.24 bioconda::samtools=1.15' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'ttps://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(reads)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta), path('*.{bam,bam.bai}'),
    path(depths), path(flagstat), path(idxstats), path(stats), emit: alignment

    script:
    def software = getSoftwareName(task.process)
    bam = "${sample_name}.Segment_${segment}.${id}.bam"
    depths = "${sample_name}.Segment_${segment}.${id}.depths.tsv"
    flagstat = "${sample_name}.Segment_${segment}.${id}.flagstat"
    idxstats = "${sample_name}.Segment_${segment}.${id}.idxstats"
    stats = "${sample_name}.Segment_${segment}.${id}.stats"
    """
    minimap2 \\
        -ax map-ont \\
        -t${task.cpus} \\
        $fasta \\
        $reads \\
    | samtools sort -@${task.cpus} \\
    > $bam
    samtools index $bam
    samtools stats $bam > $stats
    samtools flagstat $bam > $flagstat
    samtools idxstats $bam > $idxstats
    samtools depth -a -d 0 $bam > $depths
    """
}
