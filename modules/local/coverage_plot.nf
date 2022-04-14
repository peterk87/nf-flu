process COVERAGE_PLOT{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_low'
    publishDir "${params.outdir}/coverage_plots/$sample_name",
         mode: params.publish_dir_mode

    conda (params.enable_conda ? 'conda-forge::python=3.9 conda-forge::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(filt_vcf)
    val (low_coverage)

    output:
    path('*.pdf'), emit: coverage_plot

    script:
    plot_filename = "coverage_plot-${sample_name}-Segment_${segment}-${id}.pdf"
    log_scale_plot_filename = "coverage_plot-${sample_name}-Segment_${segment}-${id}-log_scale.pdf"
    """
    plot_coverage.py -d $depths -v $filt_vcf -o $plot_filename --low-coverage $low_coverage --sample-name $sample_name --segment $segment
    plot_coverage.py -d $depths -v $filt_vcf -o $log_scale_plot_filename --low-coverage $low_coverage --sample-name $sample_name --segment $segment --log-scale-y
    """
}
