// Import generic module functions
include { getSoftwareName } from './functions'

process CLAIR3{
    tag "$sample_name - Segment:$segment - Ref ID:$id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::clair3==0.1.10' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container 'hkubal/clair3:v0.1-r10'
    } else {
      container 'hkubal/clair3:v0.1-r10'
    }

    input:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(bam),
    path(depths)

    output:
    tuple val(sample_name), val(segment), val(id), path(fasta), path(depths), path(vcf), emit: vcf
    path (clair3_dir) , emit: output_dir
    path (clair3_log), emit: log
    path "versions.yml" , emit: versions

    script:
    def software  = getSoftwareName(task.process)
    vcf           = "${sample_name}.Segment_${segment}.${id}.clair3_variant.vcf.gz"
    clair3_log    = "${sample_name}.Segment_${segment}.${id}.clair3_variant.log"
    clair3_dir    = "${sample_name}.Segment_${segment}.${id}.clair3_variant"
    model_suffix  = "models/${params.clair3_variant_model}"
    """
    clair3_path=\$(which run_clair3.sh | sed 's/run_clair3.sh//g')
    if [ ${params.enable_conda} = true ] ; then
        model_path=\$clair3_path"$model_suffix"
    else
        model_path=/opt/models/${params.clair3_variant_model}
    fi
    echo "run_clair3.sh path: \$clair3_path"
    echo "models path: \$model_path"
    samtools faidx $fasta
    run_clair3.sh \\
        --bam_fn=${bam[0]} \\
        --ref_fn=$fasta \\
        --model_path="\$model_path"\\
        --threads=${task.cpus} \\
        --platform="ont" \\
        --output=${clair3_dir} \\
        --haploid_sensitive \\
        --enable_long_indel \\
        --fast_mode \\
        --include_all_ctgs
    cp ${clair3_dir}/merge_output.vcf.gz ${vcf}
    ln -s .command.log ${clair3_log}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(cat ${clair3_dir}/run_clair3.log 2>&1 | head -n1 | sed 's/^.*VERSION: //; s/ .*\$//')
    END_VERSIONS
    """
}
