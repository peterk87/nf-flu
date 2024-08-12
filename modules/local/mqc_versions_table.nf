process MQC_VERSIONS_TABLE {
    executor 'local'
    memory 1.MB

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mqc_versions_table \\
      --versions-yaml $versions \\
      --nextflow-version ${workflow.nextflow.version} \\
      --workflow-name ${workflow.manifest.name} \\
      --workflow-version ${workflow.manifest.version}
    """
}
