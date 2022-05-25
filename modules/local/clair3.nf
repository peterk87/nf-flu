// Import generic module functions
include { getSoftwareName; fluPrefix } from './functions'

process CLAIR3 {
  tag "$sample|$segment|$ref_id"
  label 'process_low'

  conda (params.enable_conda ? 'bioconda::clair3==0.1.10' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'hkubal/clair3:v0.1-r10'
  } else {
    container 'hkubal/clair3:v0.1-r10'
  }

  input:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(bam)

  output:
  tuple val(sample), val(segment), val(ref_id), path(ref_fasta), path(vcf), emit: vcf
  path (clair3_dir), emit: output_dir
  path (clair3_log), emit: log
  path "versions.yml" , emit: versions

  script:
  def software = getSoftwareName(task.process)
  def prefix   = fluPrefix(sample, segment, ref_id)
  vcf          = "${prefix}.clair3.vcf.gz"
  clair3_dir   = "${prefix}.clair3/"
  clair3_log   = "${clair3_dir}run_clair3.log"
  model_suffix = "models/${params.clair3_variant_model}"
  """
  CLAIR_BIN_DIR=\$(dirname \$(which run_clair3.sh))
  if [ ${params.enable_conda} = true ] ; then
      MODEL_PATH="\$CLAIR_BIN_DIR/${model_suffix}"
  else
      MODEL_PATH="/opt/models/${params.clair3_variant_model}"
  fi

  samtools faidx $ref_fasta

  run_clair3.sh \\
      --bam_fn=${bam[0]} \\
      --ref_fn=$ref_fasta \\
      --model_path="\$MODEL_PATH"\\
      --threads=${task.cpus} \\
      --platform="ont" \\
      --output=${clair3_dir} \\
      --haploid_sensitive \\
      --enable_long_indel \\
      --fast_mode \\
      --include_all_ctgs

  ln -s ${clair3_dir}/merge_output.vcf.gz ${vcf}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      clair3: \$(head -n1 ${clair3_dir}/run_clair3.log | sed 's/^.*CLAIR3 VERSION: v//; s/ .*\$//')
      samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
  END_VERSIONS
  """
}
