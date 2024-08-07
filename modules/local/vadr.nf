process VADR {
  tag "$sample"
  label 'process_low'

  conda 'pkru22::vadr=1.6.3'
  container 'staphb/vadr:1.6.3-hav-flu2'

  input:
  tuple val(sample), path(fasta)

  output:
  tuple val(sample), path("${prefix}/*.vadr.pass.tbl"), optional: true, emit: feature_table
  tuple val(sample), path("${prefix}/*.vadr.pass.fa"), optional: true, emit: pass_fasta
  tuple val(sample), path("${prefix}/"), emit: vadr_outdir
  path "versions.yml", emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${sample}"
  """
  v-annotate.pl \\
    $args \\
    $fasta \\
    ${prefix}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      vadr: \$(v-annotate.pl -h | perl -ne 'print "\$1\\n" if /^# VADR (\\d+\\.\\d+\\.\\d+)/')
  END_VERSIONS
  """
}

process VADR_SUMMARIZE_ISSUES {
  executor 'local'
  memory 100.MB

  input:
  path(vadr_output, stageAs: "input*/*")

  output:
  path('vadr-annotation-issues.txt'), emit: issues
  path('vadr-annotation-failed-sequences.txt'), emit: failed

  script:
  """
  cat input*/**/*.alt.list | awk 'NR == 1 || \$0 !~ /^#/' > vadr-annotation-issues.txt
  cat input*/**/*.fail.list > vadr-annotation-failed-sequences.txt
  """
}
