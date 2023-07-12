include { fluPrefix } from './functions'

process SNPEFF_ANN {
    tag "$sample|$segment|$ref_id"
    label 'process_low'

    conda 'bioconda::snpeff=5.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/snpeff:5.0--0'
    } else {
        container 'quay.io/biocontainers/snpeff:5.0--0'
    }

    input:
    tuple val(sample),
          val(segment),
          path(vcf),
          val(ref_id),
          path(fasta),
          path(db),
          path(config)

    output:
    tuple val(sample),
          val(segment),
          val(ref_id),
          path(fasta),
          path("*.vcf"), emit: vcf
    tuple val(sample),
          val(segment),
          val(ref_id),
          path(fasta),
          path("*.csv"), emit: csv
    tuple val(sample),
          val(segment),
          val(ref_id),
          path(fasta),
          path("*.genes.txt"), emit: txt
    tuple val(sample),
          val(segment),
          val(ref_id),
          path(fasta),
          path("*.html"), emit: html
    path('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix    = fluPrefix(sample, segment, ref_id)
    def avail_mem = 4
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    snpEff \\
        -Xmx${avail_mem}g \\
        ${ref_id} \\
        -config $config \\
        -dataDir $db \\
        $args \\
        $vcf \\
        -csvStats ${prefix}.snpeff.csv \\
        > ${prefix}.snpeff.vcf
    mv snpEff_summary.html ${prefix}.snpeff.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(snpEff -version 2>&1 | sed 's/^.*SnpEff //; s/ .*\$//')
    END_VERSIONS
    """
}
