include { getSoftwareName; fluPrefix } from './functions'

process MOSDEPTH_GENOME {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::mosdepth=0.3.8'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0'
  } else {
    container 'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(fasta), path(bam_bai)

  output:
  tuple val(sample), val(segment), val(ref_id), path("*.per-base.bed.gz"), emit: bedgz
  path "*.global.dist.txt", emit: mqc
  path "*.{txt,gz,csi,tsv}"
  path  "versions.yml"                          , emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix = fluPrefix(sample, segment, ref_id)
  """
  mosdepth \\
      --fast-mode \\
      $prefix \\
      ${bam_bai[0]}
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
  END_VERSIONS
  """
}
