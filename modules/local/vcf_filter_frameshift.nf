process VCF_FILTER_FRAMESHIFT{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_low'
    publishDir "${params.outdir}/variants/$sample_name",
        pattern: "*.vcf",
           mode: params.publish_dir_mode

    conda (params.enable_conda ? 'conda-forge::python=3.9 conda-forge::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(bcf_filt_vcf)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(depths), path(filt_vcf), emit: vcf

    script:
    filt_vcf = "${sample_name}.Segment_${segment}.${id}.filt.vcf"
    """
    vcf_filter_frameshift.py $bcf_filt_vcf $filt_vcf
    """
}
