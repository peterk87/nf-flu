include { getSoftwareName } from './functions'

process MOSDEPTH_GENOME {
  tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
  label 'process_medium'
  publishDir "${params.outdir}/mosdepth/$sample_name",
         mode: params.publish_dir_mode

  conda (params.enable_conda ? 'bioconda::mosdepth=0.3.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/mosdepth:0.3.1--ha7ba039_0"
  } else {
      container "quay.io/biocontainers/mosdepth:0.3.1--ha7ba039_0"
  }
  input:
  tuple val(sample_name), val(segment), val(id), path(fasta), path(bam),
    path(depths), path(flagstat), path(idxstats), path(stats)

  output:
  tuple val(sample_name), val(segment), val(id), path("*.per-base.bed.gz"), emit: bedgz
  path "*.global.dist.txt", emit: mqc
  path "*.{txt,gz,csi,tsv}"
  path '*.version.txt'                 , emit: version

  script:
  def software = getSoftwareName(task.process)
  def prefix = "${sample_name}.Segment_${segment}.${id}"
  """
  mosdepth \\
      --fast-mode \\
      $prefix \\
      ${bam[0]}
  echo \$(mosdepth --version 2>&1) | sed 's/^.*mosdepth //' > ${software}.version.txt
  """
}
