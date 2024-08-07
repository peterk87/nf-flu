// Import generic module functions
include { getSoftwareName; fluPrefix } from './functions';

process BCF_CONSENSUS {
  tag "$sample|$segment|$ref_id"
  label 'process_medium'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id) , path(fasta), path(vcf), path(mosdepth_per_base)
  val(low_coverage)

  output:
  tuple val(sample), path(consensus), emit: fasta
  path "versions.yml" , emit: versions

  script:
  def prefix = fluPrefix(sample, segment, ref_id)
  consensus    = "${prefix}.bcftools.consensus.fasta"
  sequenceID   = "${sample}_${segment}"
  """
  bgzip -c $vcf > ${vcf}.gz
  tabix ${vcf}.gz

  # get low coverage depth mask BED file by filtering for regions with less than ${low_coverage}X
  zcat $mosdepth_per_base | awk '\$4<${low_coverage}' > low_cov.bed

  bcftools consensus \\
    -f $fasta \\
    -m low_cov.bed \\
    ${vcf}.gz > $consensus

  sed -i -E "s/^>(.*)/>$sequenceID/g" $consensus

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}

process BCF_FILTER {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(fasta), path(vcf)
  val(allele_fraction)

  output:
  tuple val(sample), val(segment), val(ref_id), path(fasta), path(bcftools_filt_vcf), emit: vcf
  path "versions.yml" , emit: versions

  script:
  def exclude = (params.variant_caller == 'medaka') ? "AF < $allele_fraction" : "FILTER='RefCall' | AF < $allele_fraction"
  def prefix = fluPrefix(sample, segment, ref_id)
  bcftools_filt_vcf = "${prefix}.bcftools_filt.vcf"
  """
  bcftools norm \\
    --check-ref w \\
    -Ov \\
    -m- \\
    -f $fasta \\
    $vcf \\
    > norm.vcf

  # filter for major alleles
  bcftools filter \\
    -e "$exclude" \\
    norm.vcf \\
    -Ov \\
    -o $bcftools_filt_vcf

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}

process BCFTOOLS_STATS {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda 'bioconda::bcftools=1.20 conda-forge::gsl=2.7'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  } else {
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
  }
  input:
  tuple val(sample), val(segment), val(ref_id), path(fasta), path(vcf)

  output:
  path("*.bcftools_stats.txt"), emit: stats
  path "versions.yml" , emit: versions

  script:
  def prefix = fluPrefix(sample, segment, ref_id)
  """
  bcftools stats -F $fasta $vcf > ${prefix}.bcftools_stats.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}
