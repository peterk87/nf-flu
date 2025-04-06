process GENOFLU {
  tag "$sample"
  label 'process_low'

  conda 'bioconda::genoflu=1.06'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/genoflu:1.06--hdfd78af_0'
  } else {
    container 'quay.io/biocontainers/genoflu:1.06--hdfd78af_0'
  }

  input:
  tuple val(sample), path(consensus_fasta)

  output:
  tuple val(sample), path("${sample}.genoflu.tsv"), optional: true, emit: tsv
  tuple val(sample), path("${sample}.genoflu.xlsx"), optional: true, emit: xlsx
  path "versions.yml" , emit: versions

  script:
  """
  genoflu.py \\
    -f $consensus_fasta \\
    -n ${sample}

  mv ${sample}*.tsv ${sample}.genoflu.tsv
  mv ${sample}*.xlsx ${sample}.genoflu.xlsx

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      genoflu: \$(genoflu.py --version 2>&1 | sed 's/^.*version //')
  END_VERSIONS
  """
}
