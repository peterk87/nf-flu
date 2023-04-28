# CFIA-NCFAD/nf-flu

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[3.1.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.3)] - 2023-04-28

Patch release to fix issue to handle lowercase subtypes (e.g. `h1n5`) from NCBI Influenza DB.

## [[3.1.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.2)] - 2022-09-01

Patch release to fix issue when user reference sequences FASTA specified, but Channel from file is not treated as a value. Code has been reverted to use `file` Nextflow function.

## [[3.1.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.1)] - 2022-08-31

Patch release to fix issue when a user-specified sequences FASTA is provided and the FASTA is concatenated with the NCBI influenza sequences FASTA, but there is no new-line character at the end of the FASTA files. New line characters are added to the FASTA files to avoid incorrect concatenation.

## [[3.1.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.0)] - 2022-05-31

The workflow's name has been changed from `nf-iav-illumina` to `nf-flu` and the official repo for `nf-flu` will be [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu/) going forward.

* Added back `bin/fastq_dir_to_samplesheet.py` for Illumina `--input` samplesheet creation from Illumina FASTQ reads directory
* Fixed [issue #12](https://github.com/peterk87/nf-flu/issues/12). Nanopore sample sheet can specify a mix of single FASTQ files and/or directories containing FASTQ files. Different reads with the same sample name will be merged prior to analysis. FASTQs can be GZIP compressed and have the extensions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`. Updated CI tests to test for this flexible sample sheet handling.
* Switched to GitHub YAML form for bug report template from Markdown template.
* CI tests now output `results/pipeline_info/` and `.nextflow.log` as artifacts for easier debugging of issues.

## [[3.0.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.0.0)] - 2022-05-24

This is a major release adding a Nanopore influenza sequence analysis subworkflow using IRMA for initial assembly and BLAST against NCBI Influenza DB sequences and optionally, user-specified sequences to identify the top reference sequence for each segment for each sample. A standard read mapping/variant calling analysis is performed: for each sample, Nanopore reads are mapped separately against each gene segment reference sequence using Minimap2; variant calling of read alignments is performed using Clair3; depth-masked consensus sequence is generated using Bcftools. Consensus sequences are BLAST searched against NCBI Influenza (and user-specified sequences) to generate a BLAST summary report and H/N subtyping report. MultiQC is used to summarize results into an interactive HTML report.

NOTE: Read mapping/variant calling analysis has not been ported to the Illumina sequence analysis subworkflow.

## [[2.0.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/2.0.1)] - 2021-06-15

Patch release to fix issue [#5](https://github.com/CFIA-NCFAD/nf-flu/issues/5); added check that IRMA `amended_consensus/` exists before concatenation of consensus FASTA files.

## [[2.0.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/2.0.0)] - 2021-06-10

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

* **NB:** Parameter has been **updated** if both old and new parameter information is present.
* **NB:** Parameter has been **added** if just the new parameter information is present.
* **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
|:-----------|:------------|:------------|
| `blast`    | 2.9.0       | 2.10.0      |
| `irma`     | 0.6.7       | 1.2.1       |
| `python`   | 3.7.3       | 3.9.0       |

* **NB:** Dependency has been **updated** if both old and new version information is present.
* **NB:** Dependency has been **added** if just the new version information is present.
* **NB:** Dependency has been **removed** if new version information isn't present.
