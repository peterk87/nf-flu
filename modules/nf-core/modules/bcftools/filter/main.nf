// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCF_FILTER {
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/vcf/$sample_name",
          pattern: "*.filt.vcf",
          mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::bcftools=1.15.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/run https://depot.galaxyproject.org/singularity/bcftools:1.15--haf5b3da_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.15--haf5b3da_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(depths), path(vcf)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(depths), path(filt_vcf), emit: vcf

    script:
    filt_vcf = "${sample_name}.Segment_${segment}.${id}.filt.vcf"
    """
    bcftools filter \\
        -e 'AC[0] >= AC[1] || AC[1]<=2' \\
        $vcf \\
        -Ov \\
        -o $filt_vcf
    """
}
