process SUBTYPING_REPORT {
  memory { 
    // Dynamically determine how much memory is required for this task based on 
    // overall size of tabular blastn inputs. For a single input, allocate 2GB
    if (blastn_results instanceof nextflow.processor.TaskPath) {
      // single input file allocate 2GB
      "2 GB"
    } else if (blastn_results instanceof nextflow.util.BlankSeparatedList) {
      // multiple input files
      // mem reqs = half of sum of GB file sizes plus 2GB wiggle room
      mem_reqs = Math.ceil(0.5 * (blastn_results.collect { it.size() }.sum()) / (1024**3)) + 2
      "${mem_reqs} GB" 
    } else {
      // not TaskPath or BlankSeparatedList, then default 2GB for memory
      "2 GB"
    }
  }
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }
  

  input:
  path(genomeset)
  path(blastn_results)
  path(samplesheet)

  output:
  path('nf-flu-subtyping-report.xlsx'), emit: report
  path('parse_influenza_blast_results.log'), emit: log
  path "versions.yml", emit: versions

  script:
  """
  parse_influenza_blast_results.py \\
   --flu-metadata $genomeset \\
   --top ${params.max_top_blastn} \\
   --excel-report nf-flu-subtyping-report.xlsx \\
   --pident-threshold $params.pident_threshold \\
   --samplesheet $samplesheet \\
   $blastn_results

  ln -s .command.log parse_influenza_blast_results.log

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}
