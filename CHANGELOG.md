# peterk87/nf-iav-illumina

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.0.0](https://github.com/peterk87/nf-iav-illumina/releases/tag/2.0.0)] - 2021-06-10

### :warning: Major enhancements

* Samplesheet input (`--input samplesheet.csv`) replaces path to reads (`--reads "reads/*_R{1,2}_*.fastq.gz"`). Sample sheet can be tab-delimited (TSV) or CSV and must have a header line and 3 columns (sample name, FASTQ path/URL to forward reads, FASTQ path/URL to reverse reads).
* Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
* All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
* Updated minimum Nextflow version to `v21.04.0` (see [nextflow#572](https://github.com/nextflow-io/nextflow/issues/1964))
* Add IRMA params
  * `irma_module`: IRMA module (default: `FLU-utr`)
  * `keep_ref_deletions`: set consensus sequence deletion by ambiguation (i.e. replace ref seq with Ns) (default: `true`)
* Add BLAST subtyping params:
  * `pident_threshold`: % identity threshold (default: `0.85`)
  * `min_aln_length`: min alignment length (default: `50`)
* Replace Azure Pipelines CI with GitHub Actions CI
* add `nextflow_schema.json` and nf-core helper Jar file and Groovy scripts for params validation, printing help
* Use nf-core modules where possible
* Use nf-core module style for all processes
* Added usage and output docs
* Updated README

### Parameters

| Old parameter | New parameter                         |
|:--------------|:--------------------------------------|
| `--reads`     | `--input`                             |
|               | `--irma_module`                       |
|               | `--keep_ref_deletions`                |
|               | `--pident_threshold`                  |
|               | `--min_aln_length`                    |
|               | `--ncbi_influenza_fasta`              |
|               | `--ncbi_influenza_metadata`           |
|               | `--slurm_queue_size`                  |
|               | `--publish_dir_mode`                  |
|               | `--validate_params`                   |
|               | `--enable_conda`                      |
|               | `--singularity_pull_docker_container` |
|               | `--show_hidden_params`                |
|               | `--schema_ignore_params`              |

> **NB:** Parameter has been __updated__ if both old and new parameter information is present.

> **NB:** Parameter has been __added__ if just the new parameter information is present.

> **NB:** Parameter has been __removed__ if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
|:-----------|:------------|:------------|
| `blast`    | 2.9.0       | 2.10.0      |
| `irma`     | 0.6.7       | 1.2.1       |
| `python`   | 3.7.3       | 3.9.0       |

> **NB:** Dependency has been __updated__ if both old and new version information is present.

> **NB:** Dependency has been __added__ if just the new version information is present.

> **NB:** Dependency has been __removed__ if new version information isn't present.
