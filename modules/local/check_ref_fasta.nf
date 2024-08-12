process CHECK_REF_FASTA {
  tag "$fasta"
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
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
