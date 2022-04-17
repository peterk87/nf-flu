process CHECK_SAMPLE_SHEET {
  publishDir "${params.tracedir}/",
             mode: params.publish_dir_mode
  // using shiptv container since it has pandas, rich, typer installed
  conda (params.enable_conda ? 'bioconda::shiptv=0.4.0' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  path(samplesheet)

  output:
  path('samplesheet.fixed.csv')

  script:
  """
  check_sample_sheet.py $samplesheet ${params.platform} samplesheet.fixed.csv
  """
}
