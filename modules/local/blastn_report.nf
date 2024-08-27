process BLASTN_REPORT {
  tag "$meta.id"
  label 'process_low'
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }

  input:
  tuple val(meta), path(blastn_results)

  output:
  path('*.xlsx'), emit: report
  path('get_blastn_report.log'), emit: log
  path "versions.yml", emit: versions

  script:
  """
  get_blastn_report.py \\
    -b $blastn_results \\
    --min-aln-length ${params.min_aln_length} \\
    -x ${meta.id}-blastn-report.xlsx

  ln -s .command.log get_blastn_report.log

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
