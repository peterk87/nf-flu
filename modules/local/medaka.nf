// Import generic module functions
include { getSoftwareName } from './functions'

process MEDAKA {
  tag "$sample|$segment|$id"
  label 'process_low'

  conda 'bioconda::medaka=1.4.4'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0'
  } else {
    container 'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0'
  }

  input:
  tuple val(sample), val(segment), val(id), path(fasta), path(bam)

  output:
  tuple val(sample), val(segment), val(id), path(fasta), path(vcf), emit: vcf
  path(medaka_dir), emit: output_dir
  path(medaka_log), emit: log
  path "versions.yml", emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix   = "${sample}.Segment_${segment}.${id}"
  vcf          = "${prefix}.medaka_variant.vcf"
  medaka_dir   = "${prefix}.medaka_variant/"
  medaka_log   = "${prefix}.medaka_variant.log"
  """
  medaka_variant \\
    -o ${medaka_dir} \\
    -t ${task.cpus} \\
    -f $fasta \\
    -i ${bam[0]} \\
    -m ${params.medaka_variant_model} \\
    -s ${params.medaka_snp_model}
  
  medaka tools annotate \\
    ${medaka_dir}/round_1.vcf \\
    $fasta \\
    ${bam[0]} \\
    ${vcf} \\
    --dpsp
  
  ln -s .command.log $medaka_log

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      medaka: \$(medaka --version 2>&1 | sed 's/^.*medaka //')
  END_VERSIONS
  """
}
