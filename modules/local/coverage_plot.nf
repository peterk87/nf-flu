// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COVERAGE_PLOT{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/coverage_plots/$sample_name",
         mode: params.publish_dir_mode

    conda (params.enable_conda ? 'conda-forge::python=3.9 bioconda::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:0245b76162955df7c0f617a3169aca31e8ccfd27-0'
    } else {
        //Use it locally, while waiting for offically published
        container 'quay.io/biocontainers/mulled-v2-2d13e56148201455b8b00c4acfbf4cb74e19b1b1:b82bf22d349bbd2f68a05820dcb45edc8491614c-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(filt_vcf)

    output:
    path('*.pdf'), emit: coverage_plot

    script:
    plot_filename = "coverage_plot-${sample_name}-Segment_${segment}-${id}.pdf"
    log_scale_plot_filename = "coverage_plot-${sample_name}-Segment_${segment}-${id}-log_scale.pdf"
    """
    plot_coverage.py -d $depths -v $filt_vcf -o $plot_filename
    plot_coverage.py -d $depths -v $filt_vcf -o $log_scale_plot_filename
    """
}
