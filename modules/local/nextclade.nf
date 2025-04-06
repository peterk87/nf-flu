process NEXTCLADE_DATASET_GET {
  tag "$dataset_name"

  conda "bioconda::nextclade=3.12.0"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/nextclade:3.12.0--h9ee0642_0'
  } else {
    container 'quay.io/biocontainers/nextclade:3.12.0--h9ee0642_0'
  }

  input:
  val(dataset)

  output:
  tuple val(dataset), path(nextclade_dataset_dir), emit: dir
  path("versions.yml"), emit: versions

  script:
  dataset_name = dataset.name
  nextclade_dataset_dir = dataset_name.replace('/', '---')
  dataset_tag = dataset.tag ? "--tag ${dataset.tag}" : ""
  """
  nextclade dataset get \\
    --name $dataset_name \\
    $dataset_tag \\
    --output-dir $nextclade_dataset_dir

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/nextclade //')
  END_VERSIONS
  """
}

process NEXTCLADE_RUN {
  tag "$dataset_name|$sample"

  conda "bioconda::nextclade=3.12.0"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/nextclade:3.12.0--h9ee0642_0'
  } else {
    container 'quay.io/biocontainers/nextclade:3.12.0--h9ee0642_0'
  }

  input:
  tuple val(sample), path(fasta), val(dataset), path(nextclade_dataset_dir)

  output:
  tuple val(sample), val(dataset), path(nextclade_tsv), emit: tsv
  path("versions.yml"), emit: versions

  script:
  dataset_name = dataset.name
  nextclade_tsv = "${sample}.${dataset_name.replace('/', '---')}.nextclade.tsv"
  """
  nextclade run \\
    --jobs ${task.cpus} \\
    --input-dataset $nextclade_dataset_dir \\
    --output-tsv $nextclade_tsv \\
    $fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/nextclade //')
  END_VERSIONS
  """
}

process AGG_NEXTCLADE_TSV {
  label 'process_medium'

  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }

  input:
  path(nextclade_tsv, stageAs: "nextclade_tsv/*")
  path(manifest_csv)

  output:
  path("nextclade.tsv"), emit: tsv
  path("versions.yml"), emit: versions

  script:
  """
  agg_nextclade_tsv.py nextclade_tsv/ $manifest_csv nextclade.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      agg_nextclade_tsv.py: \$(agg_nextclade_tsv.py --version)
  END_VERSIONS
  """
}
