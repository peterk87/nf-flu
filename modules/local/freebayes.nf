// Import generic module functions
include { getSoftwareName; fluPrefix } from './functions'

process FREEBAYES {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::freebayes==1.2.0 bioconda::samtools=1.20'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-dc422f74572321627996aa5f1ae09bb26455f9de:c15b11f21039de55fb91d4e66ff0105bf3c73d7c-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-dc422f74572321627996aa5f1ae09bb26455f9de:c15b11f21039de55fb91d4e66ff0105bf3c73d7c-0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(bam)

  output:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(vcf), emit: vcf
  path (freebayes_dir), emit: output_dir
  path "versions.yml", emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix   = fluPrefix(sample, segment, ref_id)
  vcf          = "${prefix}.freebayes.vcf"
  freebayes_dir   = "${prefix}.freebayes"
  """
  mkdir -p ${freebayes_dir}

  samtools faidx $ref_fasta

  freebayes \\
      -f $ref_fasta \\
      -b ${bam[0]} \\
      --min-alternate-fraction ${params.min_alternate_fraction} \\
      --vcf ${freebayes_dir}/${vcf}

  ln -s ${freebayes_dir}/${vcf} ${vcf}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      freebayes: \$(freebayes --version 2>&1 | sed -n 's/^version:  v\\([0-9.]*\\).*/\\1/p')
  END_VERSIONS

  echo "Content of versions.yml:"
  cat versions.yml
  """
}
