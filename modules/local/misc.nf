include { saveFiles; getSoftwareName } from './functions'

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
    path "versions.yml" , emit: versions

    script:
    """
    fastqFile=\$(ls $reads | head -1)
    if grep -q "fastq.gz" <<< "\$fastqFile"; then
        cat $reads/*.fastq.gz > ${meta.id}.fastq.gz
    else
        cat $reads/*.fastq | pigz -ck > ${meta.id}.fastq.gz
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}

process GUNZIP {
    tag "$archive"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path archive

    output:
    path "$gunzip",       emit: gunzip
    path "versions.yml" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    gunzip       = archive.toString() - '.gz'
    """
    gunzip -f $archive
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

