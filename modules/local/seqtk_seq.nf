// Import generic module functions
include { getSoftwareName } from './functions'

process SEQTK_SEQ{
    tag "$sample_name - Segment:$segment - Ref ID:$id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    tuple val(sample_name), val(segment), val(id), path(reads)
    path (db)

    output:
    tuple val(sample_name), val(segment), val(id), path('*.fasta'), path(reads), emit: sample_info
    path "versions.yml"               , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = "${sample_name}"
    """
    echo $id > ${prefix}.Segment_${segment}.${id}.list
    seqtk subseq $db ${prefix}.Segment_${segment}.${id}.list > ${prefix}.Segment_${segment}.${id}.reference.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}