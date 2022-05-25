include { fluPrefix } from './functions'

process COVERAGE_PLOT{
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda (params.enable_conda ? 'python=3.9 conda-forge::typer=0.3.2 conda-forge::rich=10.6.0 conda-forge::seaborn=0.11.0 conda-forge::pandas=1.3.0 bioconda::bcbio-gff=0.6.6 bioconda::dna_features_viewer=3.0.3' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-596f42d854e849eb773ecd1b48f2b698c2d09c9f:400d0a2593841aa0bfa3402fe85debd55a29cf37-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-596f42d854e849eb773ecd1b48f2b698c2d09c9f:400d0a2593841aa0bfa3402fe85debd55a29cf37-0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(vcf), path(perbase_bed_gz)
  val(low_coverage)

  output:
  path('*.pdf'), emit: coverage_plot
  path "versions.yml", emit: versions

  script:
  def prefix = fluPrefix(sample, segment, ref_id)
  plot_filename = "coverage_plot-${prefix}.pdf"
  log_scale_plot_filename = "coverage_plot-${prefix}-log_scale.pdf"
  """
  # linear scale coverage plot
  plot_coverage.py \\
    -d $perbase_bed_gz \\
    -v $vcf \\
    -o $plot_filename \\
    --low-coverage $low_coverage \\
    --sample-name $sample \\
    --segment $segment \\
    --ref-id $ref_id

  # y-axis log scale coverage plot
  plot_coverage.py \\
    -d $perbase_bed_gz \\
    -v $vcf \\
    -o $log_scale_plot_filename \\
    --low-coverage $low_coverage \\
    --sample-name $sample \\
    --segment $segment \\
    --ref-id $ref_id \\
    --log-scale-y

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
