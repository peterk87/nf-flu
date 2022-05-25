process VCF_FILTER_FRAMESHIFT {
  tag "$sample_name|$segment|$id"
  label 'process_low'

  conda (params.enable_conda ? 'conda-forge::python=3.9 conda-forge::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0'
  }

  input:
  tuple val(sample_name), val(segment), val(ref_id), path(ref_fasta), path(vcf)

  output:
  tuple val(sample_name), val(segment), val(ref_id), path(ref_fasta), path(no_frameshifts_vcf), emit: vcf
  path "versions.yml", emit: versions

  script:
  no_frameshifts_vcf = "${sample_name}.Segment_${segment}.${ref_id}.no_frameshifts.vcf"
  """
  vcf_filter_frameshift.py $vcf $no_frameshifts_vcf

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
