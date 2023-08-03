process BLASTN_REPORT {
  tag "$meta.id"
  label 'process_low'
  conda 'conda-forge::python=3.9 conda-forge::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
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
