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

  awk '/^>/ {print; next} {gsub(/[RYSWKMBDHVryswkmbdhv]/, "N"); print}' $fasta > ${fasta}.no_ambiguous.fasta


  bcftools consensus \\
    -f ${fasta}.no_ambiguous.fasta \\
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
  val(major_allele_fraction)
  val(minor_allele_fraction)

  output:
  tuple val(sample), val(segment), val(ref_id), path(fasta), path(bcftools_filt_vcf), emit: vcf
  path "versions.yml" , emit: versions

  script:
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

  bcftools +fill-tags \\
    norm.vcf \\
    -Ov \\
    -o filled.vcf \\
    -- -t all

  bcftools +setGT \\
    filled.vcf \\
    -Ov \\
    -o setGT.major.vcf \\
    -- -t q -n 'c:1/1' -i 'FMT/VAF >= ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.major.vcf \\
    -Ov \\
    -o setGT.minor.vcf \\
    -- -t q  -n 'c:0/1' -i 'FMT/VAF >= ${minor_allele_fraction} && FMT/VAF < ${major_allele_fraction}'

  bcftools +setGT \\
    setGT.minor.vcf \\
    -Ov \\
    -o $bcftools_filt_vcf \\
    -- -t q -n 'c:0/0' -i 'FMT/VAF < ${minor_allele_fraction}'

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
