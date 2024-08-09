process PRE_TABLE2ASN {
  tag "$sample"
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }

  input:
  tuple val(sample), path(feature_table, stageAs: "input*/*"), path(fasta, stageAs: "input*/*")

  output:
  tuple val(sample), path("${sample}.tbl"), path("${sample}.fa"), path("${sample}.namesub.txt"), emit: table2asn_input
  path "versions.yml", emit: versions

  script:
  """
  sub_seqids_for_table2asn.py -t $feature_table -f $fasta -o ./ -p $sample

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}

process TABLE2ASN {
  tag "$sample"
  conda 'bioconda::table2asn=1.28.943'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/table2asn:1.28.943--h48fe88c_0'
  } else {
    container 'quay.io/biocontainers/table2asn:1.28.943--h48fe88c_0'
  }

  input:
  tuple val(sample), path(feature_table), path(fasta, stageAs: 'input*/*'), path(namesub)

  output:
  tuple val(sample), path("${prefix}.gbf"), path(namesub), optional: true, emit: genbank
  path "versions.yml", emit: versions

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${sample}"
  """
  # ensure that Genbank output file has the right name
  ln -s $fasta ${prefix}.fa
  # suppress non-zero exit codes due to warnings from table2asn
  table2asn -indir ./ -outdir ./ -V b -c - -i ${prefix}.fa || true

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      table2asn: \$(table2asn -version | sed 's/table2asn: //')
  END_VERSIONS
  """
}

process POST_TABLE2ASN {
  tag "$sample"
  conda 'bioconda::gfflu=0.0.2'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/gfflu:0.0.2--pyhdfd78af_0'
  } else {
    container 'quay.io/biocontainers/gfflu:0.0.2--pyhdfd78af_0'
  }

  input:
  tuple val(sample), path(genbank, stageAs: 'input*/*'), path(namesub, stageAs: 'input*/*')

  output:
  tuple val(sample), path("${sample}.gbk"), emit: genbank
  tuple val(sample), path("${sample}.gff"), emit: gff
  tuple val(sample), path("${sample}.faa"), emit: cds_aa_fasta
  tuple val(sample), path("${sample}.ffn"), emit: cds_nt_fasta
  path "versions.yml", emit: versions

  script:
  """
  post_table2asn.py $genbank $namesub -p $sample -o ./

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
