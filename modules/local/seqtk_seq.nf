// Import generic module functions
include { fluPrefix } from './functions'

process SEQTK_SEQ{
  tag "$sample|$segment|$ref_id"
  // use default process resources

  conda "bioconda::seqtk=1.4"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2'
  } else {
    container 'quay.io/biocontainers/seqtk:1.4--he4a0461_2'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(reads)
  path (fasta)

  output:
  tuple val(sample), val(segment), val(ref_id), path('*.fasta'), path(reads), emit: sample_info
  path "versions.yml", emit: versions

  script:
  def prefix   = fluPrefix(sample, segment, ref_id)
  """
  echo $ref_id | seqtk subseq $fasta - > ${prefix}.reference.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
  END_VERSIONS
  """
}
