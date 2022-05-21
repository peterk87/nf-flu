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
    fastqFile=\$(ls $reads | grep -E ".fastq|fq" | head -1)
    echo "\$fastqFile"
    if grep -qE ".fastq.gz" <<< "\$fastqFile"; then
        cat $reads/*.fastq.gz > ${meta.id}.fastq.gz
    elif grep -qE ".fastq" <<< "\$fastqFile"; then
        cat $reads/*.fastq | pigz -ck > ${meta.id}.fastq.gz
    elif grep -qE ".fq.gz" <<< "\$fastqFile"; then
       cat $reads/*.fq.gz > ${meta.id}.fastq.gz
    else
       cat $reads/*.fq | pigz -ck > ${meta.id}.fastq.gz
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}

process CAT_DB {
    tag "$fasta1 - $fasta2"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pigz=2.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-2b04072095278721dc9a5772e61e406f399b6030:7c7abf911e92d7fb831611ffb965f3cf7fe2c01d-0"
    }

    input:
    path(fasta1)
    path(fasta2)

    output:
    path("*.fasta"), emit: fasta

    script:
    """
    cat $fasta1 $fasta2 > influenza_db.fasta
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
    path "*.fna",       emit: gunzip
    path "versions.yml" , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    zcat $archive | sed -E 's/^>gi\\|[0-9]+\\|gb\\|(\\w+)\\|(.*)/>\\1 \\2/' > influenza.fna
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zcat: \$(echo \$(zcat --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

process CAT_CONSENSUS {
  tag "$sample_name"
  conda (params.enable_conda ? 'bioconda::shiptv=0.4.0' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  tuple val(sample_name), path(consensus)

  output:
  tuple val(sample_name), path('*.consensus.blastn.fasta'), emit: fasta
  path('*.consensus.fasta'), emit: consensus_fasta
  path "versions.yml" , emit: versions

  script:
  """
  cat_consensus_sequences.py \\
  --sample-name $sample_name \\
  --output1-fasta ${sample_name}.consensus.fasta \\
  --output2-fasta ${sample_name}.consensus.blastn.fasta \\
  $consensus
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}



