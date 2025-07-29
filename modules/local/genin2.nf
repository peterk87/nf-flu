process GENIN2 {
  label 'process_low'

  conda 'bioconda::genin2=2.1.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/genin2:2.1.3--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/genin2:2.1.3--pyhdfd78af_0'
  }

  input:
  path(fasta)

  output:
  path('genin2.tsv'), emit: tsv
  path('versions.yml'), emit: versions

  script:
  """
  genin2 -o genin2.tsv $fasta

  #> genin --version
  #  genin2, version 2.1.3, by Alessandro Sartori (asartori@izsvenezie.it)
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      genin2: \$(genin2 --version 2>&1 | sed 's/^.*genin2, version //' | sed 's/,.*//')
  """
}