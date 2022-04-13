process BCF_CONSENSUS {
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/consensus/bcf_consensus/$sample_name",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::bcftools=1.15.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/run https://depot.galaxyproject.org/singularity/bcftools:1.15--haf5b3da_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.15--haf5b3da_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id) , path(fasta), path(vcf), path(mosdepth_per_base)
    val(low_coverage)

    output:
    tuple val(sample_name), val(segment), val(id), path(consensus), emit: vcf

    script:
    consensus = "${sample_name}.Segment_${segment}.${id}.consensus.fasta"
    sequenceID = "${sample_name}.Segment${segment}.${id}"
    """
    bgzip -c $vcf > ${vcf}.gz
    tabix ${vcf}.gz
    zcat $mosdepth_per_base | awk '\$4<${low_coverage}' > low_cov.bed
    bcftools consensus \\
        -f $fasta \\
        -m low_cov.bed \\
        ${vcf}.gz > $consensus
    sed -i -E "s/^>(.*)/>$sequenceID/g" $consensus
  """
}

process BCF_FILTER {
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.15.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/run https://depot.galaxyproject.org/singularity/bcftools:1.15--haf5b3da_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.15--haf5b3da_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(vcf)
    val (allele_fraction)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta),
    path(depths), path(filt_vcf), emit: vcf

    script:
    filt_vcf = "${sample_name}.Segment_${segment}.${id}.filt.vcf"
    """
    bcftools norm \\
        -Ov \\
        -m- \\
        -f $fasta \\
        $vcf \\
        > norm.vcf
    # filter for major alleles
    bcftools filter \\
        -e 'AF < $allele_fraction' \\
        norm.vcf \\
        -Ov \\
        -o $filt_vcf
    """
}