// Import generic module functions
include { getSoftwareName } from './functions'

process MEDAKA{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/variants/$sample_name",
        pattern: "*.{vcf,log}",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::medaka=1.4.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container 'https://depot.galaxyproject.org/singularity/medaka:1.4.3--py38h130def0_0'
    } else {
      container 'quay.io/biocontainers/medaka:1.4.4--py38h130def0_0'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(bam),
    path(depths)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(vcf), emit: vcf
    path('*.medaka.log'), emit: log
    path '*.version.txt'                 , emit: version

    script:
    def software  = getSoftwareName(task.process)
    vcf           = "${sample_name}.Segment_${segment}.${id}.medaka.vcf"
    medaka_log    = "${sample_name}.Segment_${segment}.${id}.medaka.log"
    """
    medaka_variant \\
        -d \\
        -o medaka_variant \\
        -t ${task.cpus} \\
        -f $fasta \\
        -i ${bam[0]} \\
        -m ${params.medaka_variant_model} \\
        -s ${params.medaka_snp_model}
    medaka tools annotate \\
        medaka_variant/round_1.vcf \\
        $fasta \\
        ${bam[0]} \\
        ${vcf} \\
        --dpsp
    ln -s .command.log $medaka_log
    echo \$(medaka --version 2>&1) | sed 's/^.*medaka //' > ${software}.version.txt
    """
}
