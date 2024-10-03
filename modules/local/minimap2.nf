include { fluPrefix } from './functions'

process MINIMAP2 {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::minimap2=2.28 bioconda::samtools=1.20'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-058de387f9917a7a63953f496cdd203bca83b790:86215829f86df9201683956877a19d025261ff66-0'
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

  // Determine the mapping option based on the platform
  def map_option = params.platform == 'nanopore' ? 'map-ont' : 'sr'
  samtools_view_flag = params.output_unmapped_reads ? "" : "-F 4"
  """
  minimap2 \\
    -ax $map_option \\
    -t${task.cpus} \\
    $ref_fasta \\
    $reads \\
    | samtools sort -@${task.cpus} \\
    | samtools view -b $samtools_view_flag \\
    > $bam

  samtools index $bam

  samtools stats $bam > $stats
  samtools flagstat $bam > $flagstat
  samtools idxstats $bam > $idxstats

  ln -s .command.log $minimap2_log

  ${params.platform == 'illumina' ? "samtools faidx $ref_fasta" : ""}
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      minimap2: \$(minimap2 --version 2>&1)
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
