process CHECK_REF_FASTA {
  tag "$fasta"
  conda 'bioconda::shiptv=0.4.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  path(fasta)

  output:
  path('ref_fasta.fixed.fasta'), emit: fasta
  path "versions.yml", emit: versions

  script:
  """
  ref_fasta_check.py $fasta ref_fasta.fixed.fasta
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
