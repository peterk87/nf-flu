process GENOFLU {
  label 'process_low'

  conda 'bioconda::genoflu=1.05'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/genoflu:1.05--hdfd78af_0'
  } else {
    container 'quay.io/biocontainers/genoflu:1.05--hdfd78af_0'
  }

  input:
  tuple val(id), path(consensus_fasta)

  output:
  tuple val(id), path("${id}*tsv"), optional: true, emit: tsv
  tuple val(id), path("${id}*xlsx"), optional: true, emit: xlsx
  path "versions.yml" , emit: versions

  script:
  """
  genoflu.py \\
    -f $consensus_fasta \\
    -n ${id}

  mv ${id}*.tsv ${id}_genoflu.tsv
  mv ${id}*.xlsx ${id}_genoflu.xlsx

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      genoflu: \$(genoflu.py --version 2>&1 | sed 's/^.*genoflu\.py: version //;s/, .*//')
  END_VERSIONS
  """
}
