process PREP_FLUMUT_FASTA {
  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }

  input:
  path(fastas)

  output:
  path('seqs-for-flumut.fasta'), emit: fasta

  script:
  """
  awk 'BEGIN {FS=OFS=""}
     /^>/ {
         gsub(/_segment7_M/, "_MP", \$0);  # Replace "_segment7_M" with "_MP"
         gsub(/segment[0-9]+_/, "", \$0);  # Remove "segment1-9_"
     }
     {print}' $fastas > seqs-for-flumut.fasta
  """
}

process FLUMUT {
  label 'process_low'

  conda 'bioconda::flumut=0.6.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/flumut:0.6.3--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/flumut:0.6.3--pyhdfd78af_0'
  }

  input:
  path(fasta)

  output:
  path('flumut.xlsm'), emit: xlsm
  path('flumut-markers.tsv'), emit: markers
  path('flumut-mutations.tsv'), emit: mutations
  path "versions.yml" , emit: versions

  script:
  """
  flumut \\
    -x flumut.xlsm \\
    -m flumut-markers.tsv \\
    -M flumut-mutations.tsv \\
    $fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      flumut: \$(flumut --version 2>&1 | sed 's/^.*flumut, version //;s/, .*//')
  """
}
