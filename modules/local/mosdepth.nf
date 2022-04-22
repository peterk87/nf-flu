include { getSoftwareName } from './functions'

process MOSDEPTH_GENOME {
  tag "$sample_name - Segment:$segment - Ref ID:$id"
  label 'process_low'

  conda (params.enable_conda ? 'bioconda::mosdepth=0.3.3' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/mosdepth:0.3.3--h37c5b7d_2"
  } else {
      container "quay.io/biocontainers/mosdepth:0.3.3--h01d7912_0"
  }
  input:
  tuple val(sample_name), val(segment), val(id), path(fasta), path(bam),
    path(depths), path(flagstat), path(idxstats), path(stats)

  output:
  tuple val(sample_name), val(segment), val(id), path("*.per-base.bed.gz"), emit: bedgz
  path "*.global.dist.txt", emit: mqc
  path "*.{txt,gz,csi,tsv}"
  path  "versions.yml"                          , emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix = "${sample_name}.Segment_${segment}.${id}"
  """
  mosdepth \\
      --fast-mode \\
      $prefix \\
      ${bam[0]}
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
  END_VERSIONS
  """
}
