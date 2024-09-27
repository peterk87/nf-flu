// Import generic module functions
include { getSoftwareName; fluPrefix } from './functions'

process FREEBAYES {
    tag "$sample|$segment|$ref_id"
    label 'process_low'

    // Container settings for different environments
    conda 'bioconda::freebayes==1.3.8'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/freebayes:1.3.8--h6a68c12_2'
    } else {
        container 'quay.io/biocontainers/freebayes:1.3.8--h6a68c12_1'
    }

    input:
    tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(bam)


    output:
    tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(vcf), emit: vcf
    path (freebayes_dir), emit: output_dir
    path "versions.yml", emit: versions

    script:
    def prefix   = fluPrefix(sample, segment, ref_id)
    vcf          = "${prefix}.freebayes.vcf"
    freebayes_dir   = "${prefix}.freebayes"
    """
    mkdir -p ${freebayes_dir}

    freebayes \\
        -f $ref_fasta \\
        -b ${bam[0]} \\
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