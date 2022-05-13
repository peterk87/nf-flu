# Influenza Genome Analysis Nextflow Workflow

[![CI](https://github.com/peterk87/nf-iav-illumina/actions/workflows/ci.yml/badge.svg)](https://github.com/peterk87/nf-iav-illumina/actions/workflows/ci.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-flu** is a bioinformatics analysis pipeline for assembly and H/N subtyping of Influenza virus. The pipeline supports both Illumina and Nanopore Platform.
Since Influenza is a special virus with multiple gene segments (8 segments) and there might be a reference or multiple we would want to align against, the pipeline will automatically pull top match references for each segment.
To achieve this task, the pipeline downloads Influenza database from NCBI and user could provide their own reference database. The pipline performs read mapping against each reference segment, variant calling and genome assembly.

The pipeline is implemented in [Nextflow][]

## Pipeline summary

1. Download latest [NCBI Influenza DB][] sequences and metadata (or use user-specified files)
2. Merge reads of re-sequenced samples ([`cat`](http://www.linfo.org/cat.html)) (if needed)
3. Assembly of Influenza gene segments with [IRMA][] using the built-in FLU module
4. Nucleotide [BLAST][] search against [NCBI Influenza DB][]
5. Automatically pull top match references for segments
6. H/N subtype prediction and Excel XLSX report generation based on BLAST results
7. Perform Variant calling and genome assembly for all segments.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`).
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort)_
3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run peterk87/nf-iav-illumina -profile test,<docker/singularity/podman/shifter/charliecloud/conda>
    ```

    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Run your own analysis

    * [Optional] Generate an input samplesheet from a directory containing Illumina FASTQ files (e.g. `/path/to/illumina_run/Data/Intensities/Basecalls/`) with the included Python script [`fastq_dir_to_samplesheet.py`](https://github.com/peterk87/nf-iav-illumina/blob/master/bin/fastq_dir_to_samplesheet.py) **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```bash
        python ~/.nextflow/assets/peterk87/nf-iav-illumina/bin/fastq_dir_to_samplesheet.py \
          -i /path/to/illumina_run/Data/Intensities/Basecalls/ \
          -o samplesheet.csv
        ```
    
    * Typical command for Illumina Platform

        ```bash
        nextflow run peterk87/nf-iav-illumina \
          --input samplesheet.csv \
          --platform illumina \
          --profile <docker/singularity/podman/shifter/charliecloud/conda>
        ```
    
    * Typical command for Nanopore Platform

      ```bash
      nextflow run peterk87/nf-iav-illumina \
        --input samplesheet.csv \
        --platform nanopore \
        --profile <docker/singularity/conda>
      ```

## Documentation

The nf-flu pipeline comes with:

* [usage](docs/usage.md) and
* [output](docs/output.md) documentation.

## Resources

* [NCBI Influenza FTP site](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/)
* [IRMA][] Iterative Refinement Meta-Assembler
  * [IRMA Publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6)

## Credits

* [nf-core](https://nf-co.re) project for establishing Nextflow workflow development best-practices, [nf-core tools](https://nf-co.re/tools-docs/) and [nf-core modules](https://github.com/nf-core/modules)
* [nf-core/viralrecon](https://github.com/nf-core/viralrecon) for inspiration and setting a high standard for viral sequence data analysis pipelines
* [Conda](https://docs.conda.io/projects/conda/en/latest/) and [Bioconda](https://bioconda.github.io/) project for making it easy to install, distribute and use bioinformatics software.
* [Biocontainers](https://biocontainers.pro/) for automatic creation of [Docker] and [Singularity] containers for bioinformatics software in [Bioconda]

[NCBI Influenza DB]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[Nextflow]: https://www.nextflow.io/
[Docker]: https://www.docker.com/
[Singularity]: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
