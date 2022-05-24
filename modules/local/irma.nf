// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IRMA {
  tag "$meta.id"
  label 'process_medium'
  publishDir "${params.outdir}/irma/", 
             mode: params.publish_dir_mode
  publishDir "${params.outdir}/consensus/",
             pattern: "${meta.id}.consensus.fasta",
             mode: params.publish_dir_mode

  conda (params.enable_conda ? "bioconda::irma=1.0.2" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/irma:1.0.2--pl5222hdfd78af_0'
  } else {
    container 'quay.io/biocontainers/irma:1.0.2--pl5222hdfd78af_0'
  }

  input:
  tuple val(meta), path(reads)
  
  output:
  tuple val(meta), path("${meta.id}/"), emit: irma
  tuple val(meta), path("${meta.id}.consensus.fasta"), optional: true, emit: consensus
  path "*.version.txt", emit: version

  script:
  def software = getSoftwareName(task.process)
  irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
  """
  touch irma_config.sh
  echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
  echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
  if [ ${params.keep_ref_deletions} ]; then
    echo 'DEL_TYPE="NNN"' >> irma_config.sh
    echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
  fi

  IRMA $params.irma_module $reads $meta.id
  
  if [ -d "${meta.id}/amended_consensus/" ]; then
    cat ${meta.id}/amended_consensus/*.fa > ${meta.id}.consensus.fasta
  fi

  set +e
  IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/' > ${software}.version.txt
  set -e
  """
}
