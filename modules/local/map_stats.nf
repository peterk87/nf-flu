// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MAP_STATS{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/mapping/$sample_name",
         mode: params.publish_dir_mode

    conda (params.enable_conda ? 'conda-forge::perl=5.32.1 bioconda::samtools=1.15' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-1a8bb10cd2d4d4773f9fec201ad000e3cc10174f:6928ed3d7c71a179513b4d1eba6315070be0eb6c-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-1a8bb10cd2d4d4773f9fec201ad000e3cc10174f:6928ed3d7c71a179513b4d1eba6315070be0eb6c-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(bam)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(bam), path(depths), path(flagstat), path(idxstats), emit: stat

    script:
    depths = "${sample_name}.Segment_${segment}.${id}.depths.tsv"
    flagstat = "${sample_name}.Segment_${segment}.${id}.flagstat"
    idxstats = "${sample_name}.Segment_${segment}.${id}.idxstats"
    """
    samtools flagstat $bam > $flagstat
    samtools depth -a -d 0 $bam | perl -ne 'chomp \$_; print "${sample_name}\t\$_\n"' > $depths
    samtools idxstats $bam | head -n1 | perl -ne 'chomp \$_; print "${sample_name}\t\$_\n"' > $idxstats
    """
}
