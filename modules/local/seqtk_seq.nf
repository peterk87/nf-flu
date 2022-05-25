// Import generic module functions
include { getSoftwareName } from './functions'

process SEQTK_SEQ{
  tag "$sample_name|$segment|$id"
  label 'process_low'

  conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3'
  } else {
    container 'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3'
  }

  input:
  tuple val(sample_name), val(segment), val(id), path(reads)
  path (db)

  output:
  tuple val(sample_name), val(segment), val(id), path('*.fasta'), path(reads), emit: sample_info
  path "versions.yml", emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix   = "${sample_name}.Segment_${segment}.${id}"
  """
  echo $id > ${prefix}.list
  seqtk subseq $db ${prefix}.list > ${prefix}.reference.fasta
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
  END_VERSIONS
  """
}
