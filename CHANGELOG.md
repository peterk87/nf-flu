# CFIA-NCFAD/nf-flu

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[3.3.4](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.4)] - 2023-08-18

### Fixes

* Subtyping report summary sheet "1_Subtype Predictions" shows only N subtype results

## [[3.3.3](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.3)] - 2023-08-16

This release fixes issues with subtype report generation script (`parse_influenza_blast_results.py`), primarily subtype predictions being `N/A` for samples where the top BLAST hits are user-specified sequences for the HA and NA segments.

### Fixes

* subtype prediction based off majority H/N prediction of all BLAST hits instead of just the top X matches (#40)
* the top hit for H/N can also be a user-specified sequence without subtype information
* top segment matches are now sorted by sample name, segment name and BLAST bitscore
* output concatenated Nanopore FASTQ to `${outdir}/fastq` by default (#43)
* Handle ambiguous bases in reference sequences by having Clair3 not convert those positions to N and Bcftools produce a warning instead of an error (#42)

### Changes

* subtyping report results are now ordered in the same order as the input `samplesheet.csv`, that is the order of the samples in the report is the same as the order of the samples in the `samplesheet.csv` file

## [[3.3.2](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.2)] - 2023-08-03

This patch release fixes an IBV subtype/genotype parsing issue when generating subtyping report using the new metadata format introduced in 3.3.0 ([#32](https://github.com/CFIA-NCFAD/nf-flu/issues/32)).

## [[3.3.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.1)] - 2023-08-02

### Fixes

* Conda/Mamba env creation when using `conda`/`mamba` profile (#35)

## [[3.3.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.3.0)] - 2023-07-11

This release migrates to more recently updated Influenza virus sequences since the last update for the [NCBI Influenza DB FTP data](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/) was in 2020-10-13. By default, all Orthomyxoviridae virus sequences were parsed from the daily updated NCBI Viruses [`AllNucleotide.fa`](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/) and [`AllNuclMetadata.csv.gz`](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv.gz) and uploaded to [Figshare](https://figshare.com/articles/dataset/2023-06-14_-_NCBI_Viruses_-_Orthomyxoviridae/23608782) as Zstd compressed files. nf-flu no longer uses the [influenza.fna.gz](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz) and [genomeset.dat.gz](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/genomeset.dat.gz) files for Influenza sequences and metadata, respectively.

### Fixes

* More up-to-date Influenza sequences database used by default (#24)

## [[3.2.1](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.2.1)] - 2023-07-07

### Fixes

* Empty BLAST results file parsing `NoDataError` (#27) (Thanks @MatFish for reporting this issue!)

## [[3.2.0](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.2.0)] - 2023-06-22

### Added

* Influenza B virus support (#14)
* Polars for faster parsing of BLAST results (#14)

### Fixes

* Irregular Illumina paired-end FASTQ files not producing IRMA assemblies (#20)

### Updates

* Updated README.md to include references and citations

## [[3.1.6](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.6)] - 2023-05-31

This is a patch release for a minor change to use Biocontainers Docker and Singularity images for Clair3 to avoid hitting limits on pulls from Docker Hub and since Biocontainers images are half the size of [hkubal/clair3](https://hub.docker.com/r/hkubal/clair3/) images.

Also, updated CI workflow and added issue template forms for feature request and questions.

## [[3.1.5](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.5)] - 2023-05-30

### Added

* `--use_mamba` to enable using [Mamba](https://github.com/mamba-org/mamba/) in place of Conda when using `-profile conda` for faster creation of Conda environments

### Updates

* Clair3: 0.1.10 -> 1.0.2

### Fixes

* user-specified Clair3 models not being found ([#11](https://github.com/CFIA-NCFAD/nf-flu/issues/11))
* Conda profile not enabling Conda ([#15](https://github.com/CFIA-NCFAD/nf-flu/issues/15))
* IRMA wanting too much `/tmp` space; IRMA's tmp dir will be output to the current working directory of the process job ([#13](https://github.com/CFIA-NCFAD/nf-flu/issues/13)) (Thanks @Codes1985 for reporting and solving this issue!)

## [[3.1.4](https://github.com/CFIA-NCFAD/nf-flu/releases/tag/3.1.4)] - 2023-05-17

This release addresses issue [#11](https://github.com/CFIA-NCFAD/nf-flu/issues/11) adding a new option `--clair3_user_variant_model <PATH TO CLAIR3 MODEL>` to allow user to provide a Clair3 model not included with Clair3, e.g. a [Rerio](https://github.com/nanoporetech/rerio) Clair3 model for r10 flowcells.

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
