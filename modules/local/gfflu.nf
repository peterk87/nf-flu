process GFFLU {
  conda "bioconda::gfflu=0.0.1"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container ''
  } else {
    container ''
  }

  input:
  tuple val(seqid), path(fasta, stageAs: "input*/*")

  output:
  tuple val(seqid), path(fasta), path("gfflu/${seqid}.gff"), emit: gff
  tuple val(seqid), path(fasta), path("gfflu/${seqid}.gbk"), emit: gbk
  path('gfflu/'), emit: outdir
  path('versions.yml'), emit: versions

  script:
  input_fasta = "${seqid}.fasta"
  """
  ln -s $fasta $input_fasta
  gfflu -v -o gfflu $input_fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      gfflu: \$(gfflu --version | sed 's/gfflu version //')
  END_VERSIONS
  """
}
