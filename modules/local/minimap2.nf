include { fluPrefix } from './functions'

process MINIMAP2 {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda (params.enable_conda ? 'bioconda::minimap2=2.24 bioconda::samtools=1.15' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0'
  }


  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(reads)

  output:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path('*.{bam,bam.bai}'), emit: alignment
  path '*.{flagstat,idxstats,stats}', emit: stats
  path('*.minimap2.log'), emit: log
  path "versions.yml" , emit: versions

  script:
  def prefix  = fluPrefix(sample, segment, ref_id)
  bam         = "${prefix}.bam"
  depths      = "${prefix}.depths.tsv"
  flagstat    = "${prefix}.flagstat"
  idxstats    = "${prefix}.idxstats"
  stats       = "${prefix}.stats"
  minimap2_log = "${prefix}.minimap2.log"
  """
  minimap2 \\
    -ax map-ont \\
    -t${task.cpus} \\
    $ref_fasta \\
    $reads \\
    | samtools sort -@${task.cpus} \\
    > $bam

  samtools index $bam

  samtools stats $bam > $stats
  samtools flagstat $bam > $flagstat
  samtools idxstats $bam > $idxstats

  ln -s .command.log $minimap2_log

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: \$(minimap2 --version 2>&1)
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
