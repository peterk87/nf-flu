process MINIMAP2{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::minimap2=2.24 bioconda::samtools=1.15' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(reads)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta), path('*.{bam,bam.bai}'),
    path(depths), path(flagstat), path(idxstats), path(stats), emit: alignment
    path('*.minimap2.log'), emit: log
    path '*.version.txt'                 , emit: version

    script:
    bam           = "${sample_name}.Segment_${segment}.${id}.bam"
    depths        = "${sample_name}.Segment_${segment}.${id}.depths.tsv"
    flagstat      = "${sample_name}.Segment_${segment}.${id}.flagstat"
    idxstats      = "${sample_name}.Segment_${segment}.${id}.idxstats"
    stats         = "${sample_name}.Segment_${segment}.${id}.stats"
    mapping_log   = "${sample_name}.Segment_${segment}.${id}.minimap2.log"
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
    ln -s .command.log $mapping_log
    echo \$(minimap2 --version 2>&1) > minimap2.version.txt
    """
}
