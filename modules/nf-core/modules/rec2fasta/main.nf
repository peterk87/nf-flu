// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REC2FASTA {
    tag "${record.id} - ${record.desc} - ${record.sequence.size()}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val (record)

    output:
    path("*.fasta"), emit: record

    script:
    fasta = "${record.id}.fasta"
  """
  cat > $fasta << EOF
>${record.id} ${record.desc}
${record.sequence}
EOF
  """
}
