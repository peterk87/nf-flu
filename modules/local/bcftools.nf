// Import generic module functions
include { getSoftwareName } from './functions';

process BCF_CONSENSUS {
    tag "$sample_name - Segment:$segment - Ref ID:$id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.15.1 conda-forge::gsl=2.7' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/run https://depot.galaxyproject.org/singularity/bcftools:1.15--haf5b3da_0'
    } else {
        container 'quay.io/biocontainers/bcftools:1.15--haf5b3da_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id) , path(fasta), path(vcf), path(mosdepth_per_base)
    val(low_coverage)

    output:
    tuple val(sample_name), path(consensus), emit: fasta
    path "versions.yml" , emit: versions

    script:
    consensus    = "${sample_name}.Segment_${segment}.${id}.bcftools.consensus.fasta"
    sequenceID   = "${sample_name}_${segment}"
    """
    bgzip -c $vcf > ${vcf}.gz
    tabix ${vcf}.gz
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
    tag "$sample_name - Segment:$segment - Ref ID:$id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::bcftools=1.15.1 conda-forge::gsl=2.7' : null)
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
    path(depths), path(bcftools_filt_vcf), emit: vcf
    path "versions.yml" , emit: versions

    script:
    bcftools_filt_vcf = "${sample_name}.Segment_${segment}.${id}.bcftools_filt.vcf"
    def exclude
    if (params.variant_caller == 'medaka'){
        exclude = "AF < $allele_fraction"
    }else{
        exclude = "%FILTER='RefCall' | AF < $allele_fraction"
    }
    """
    bcftools norm \\
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