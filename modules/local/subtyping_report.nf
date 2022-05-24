// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SUBTYPING_REPORT {
  publishDir "${params.outdir}/", mode: params.publish_dir_mode
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
  conda (params.enable_conda ? 'conda-forge::python=3.9 bioconda::biopython=1.78 conda-forge::openpyxl=3.0.7 conda-forge::pandas=1.2.4 conda-forge::rich=10.2.2 conda-forge::typer=0.3.2 conda-forge::xlsxwriter=1.4.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0'
    }

  input:
  path(genomeset)
  path(blastn_results)

  output:
  path('iav-subtyping-report.xlsx'), emit: report
  path('parse_influenza_blast_results.log'), emit: log

  script:
  """
  parse_influenza_blast_results.py \\
   --threads ${task.cpus} \\
   --flu-metadata $genomeset \\
   --excel-report iav-subtyping-report.xlsx \\
   --pident-threshold $params.pident_threshold \\
   $blastn_results
  ln -s .command.log parse_influenza_blast_results.log
  """
}
