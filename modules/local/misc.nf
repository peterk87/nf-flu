process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pigz=2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    """
    fastqFile=\$(ls $reads | head -1)
    if grep -q "gz" <<< "\$fastqFile"; then
        cat $reads/*.fastq.gz > ${meta.id}.fastq.gz
    else
        cat $reads/*.fastq | pigz -ck > ${meta.id}.fastq.gz
    fi
    """
}
