// modified from nf-core/modules version to only mv BLAST DB files to blast_db dir
process BLAST_MAKEBLASTDB {
  tag "$fasta"
  label 'process_low'

  conda 'bioconda::blast=2.14.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/blast:2.14.0--h7d5a4b4_1'
  } else {
    container 'quay.io/biocontainers/blast:2.14.0--h7d5a4b4_1'
  }

  input:
  path fasta

  output:
  path 'blast_db'     , emit: db
  path "versions.yml" , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
  makeblastdb \\
      -in $fasta \\
      $args
  mkdir blast_db
  mv ${fasta}.* blast_db
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
  END_VERSIONS
  """
}
