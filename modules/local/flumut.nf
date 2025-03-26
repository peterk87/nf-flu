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

process VADR2FLUMUT {
  label 'process_low'
  
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }

  input:
  path(vadr_gbks, stageAs: "vadr_gbks/*")

  output:
  path('seqs-for-flumut.fasta'), emit: fasta

  script:
  """
  vadr2flumut.py vadr_gbks/ > seqs-for-flumut.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      vadr2flumut.py: \$(subtyping_report.py --version)
  END_VERSIONS
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
