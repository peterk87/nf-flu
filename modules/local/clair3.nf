// Import generic module functions
include { getSoftwareName } from './functions'

process CLAIR3{
    tag "$sample_name - Segment:$segment - Ref ID:$id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::clair3==0.1.10--hdfd78af_0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container 'https://depot.galaxyproject.org/singularity/clair3:0.1.10--hdfd78af_0'
    } else {
      container 'quay.io/biocontainers/clair3:0.1.10--hdfd78af_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(bam),
    path(depths)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(vcf), emit: vcf
    path('*.clair3_variant.log'), emit: log
    path "versions.yml" , emit: versions

    script:
    def software  = getSoftwareName(task.process)
    vcf           = "${sample_name}.Segment_${segment}.${id}.clair3_variant.vcf.gz"
    clair3_log    = "${sample_name}.Segment_${segment}.${id}.clair3_variant.log"
    """
    samtools faidx $fasta
    run_clair3.sh \\
        --bam_fn=${bam[0]} \\
        --ref_fn=$fasta \\
        --model_path="/usr/local/bin/models/${params.clair3_variant_model}" \\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --output=clair3_variant \\
        --haploid_sensitive \\
        --enable_long_indel \\
        --fast_mode \\
        --include_all_ctgs
    cp clair3_variant/merge_output.vcf.gz ${vcf}
    ln -s .command.log $clair3_log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(cat clair3_variant/run_clair3.log 2>&1 | head -n1 | sed 's/^.*VERSION: //; s/ .*\$//')
    END_VERSIONS
    """
}
