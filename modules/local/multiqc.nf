process MULTIQC {
  label 'process_long'

  conda "bioconda::multiqc=1.23"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/multiqc:1.23--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0'
  }

  input:
  path(multiqc_custom_config)
  path('samtools/*')
  path('mosdepth/*')
  path('bcftools/*')
  path('software_versions/*')
  path(workflow_summary)

  output:
  path "*.html", emit: multiqc_report
  path "*_data"
  path "multiqc_plots"
  path "versions.yml", emit: versions

  script:
  def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
  """
  multiqc -f $custom_config .
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
  END_VERSIONS
  """
}
