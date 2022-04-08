// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BLAST_BLASTDBCMD{
    tag "$sample_name - Segment:$segment - Ref Accession ID:$id"
    label 'process_medium'
    publishDir "${params.outdir}/reference_sequences/",
        pattern: "*.fasta",
        mode: params.publish_dir_mode

    conda (params.enable_conda ? 'bioconda::blast=2.10.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/blast:2.10.1--pl526he19e7b1_3'
    } else {
        container 'quay.io/biocontainers/blast:2.10.1--pl526he19e7b1_3'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(reads)
    path (db)

    output:
    tuple val(sample_name), val(segment), val(id), path('*.fasta'), path(reads), emit: fasta
    path '*.version.txt'            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${sample_name}${options.suffix}" : "${sample_name}"
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/.ndb//'`
    blastdbcmd \\
        -db \$DB \\
        -entry $id\\
        -out ${prefix}.Segment_${segment}.${id}.reference.fasta
    echo \$(blastn -version 2>&1) | sed 's/^.*blastn: //; s/ .*\$//' > ${software}.version.txt
    """
}
