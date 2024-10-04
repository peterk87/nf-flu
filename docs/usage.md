# nf-flu: Usage

## Introduction

### `--input` samplesheet

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a tab-delimited (TSV) or comma-separated (CSV) file with 2 or 3 columns, for Nanopore or Illumina sequence data, respectively, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

#### Illumina samplesheet

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will concatenate the raw reads before performing any downstream analysis.

Example samplesheet where `SAMPLE_1` was sequenced twice with Illumina paired-end sequencing and `SAMPLE_2` was sequenced once with Illumina single-end sequencing.

```bash
sample,fastq_1,fastq_2
SAMPLE_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
SAMPLE_2,AEG588A2_S4_L003_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                 |
|-----------|-----------------------------------------------------------------------------------------------------------------------------|
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.               |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".  |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".  |

#### Generate a samplesheet from a FASTQ directory

You can generate an input samplesheet from a directory containing Illumina FASTQ files (e.g. `/path/to/illumina_run/Data/Intensities/Basecalls/`) with the included Python script [`fastq_dir_to_samplesheet.py`](bin/fastq_dir_to_samplesheet.py) **before** you run the pipeline (requires Python 3 installed locally) e.g.

```bash
~/.nextflow/assets/CFIA-NCFAD/nf-flu/bin/fastq_dir_to_samplesheet.py \
  -i /path/to/illumina_run/Data/Intensities/Basecalls/ \
  -o samplesheet.csv
```

#### Nanopore samplesheet

The samplesheet for Nanopore sequencing data analysis must be either a comma-separated (CSV) or tab-delimited (TSV) file with 2 columns and a header with column names (names don't matter).

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will concatenate the raw reads before performing any downstream analysis.

Example Nanopore samplesheet:

```bash
sample,reads
SAMPLE_1,/path/to/run1/fastq_pass/barcode01
SAMPLE_1,/path/to/ANOTHER_RUN/fastq_pass/barcode09
SAMPLE_1,/path/to/sample1.fastq.gz
SAMPLE_2,/path/to/run2/fastq_pass/barcode02
```

- **NOTE:** `SAMPLE_1` has 3 entries

| Column    | Description                                                                                                                                       |
|-----------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.                                     |
| `reads`   | Full path to FASTQ file or directory containing basecalled reads in gzipped or unzipped FASTQ file format with extension ".fastq.gz" or ".fastq". |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run CFIA-NCFAD/nf-flu --input samplesheet.csv -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull CFIA-NCFAD/nf-flu
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [CFIA-NCFAD/nf-flu releases page](https://github.com/CFIA-NCFAD/nf-flu/releases) and find the latest version number - numeric only (eg. `2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Pipeline parameters

### Input/output options

Define where the pipeline should find input data and save output data.

#### `--input`

- **Required**
- Type: string

Sample sheet with sample names and paths to reads.

#### `--platform`

- Optional
- Type: string
- Default: illumina

Specify the platform for Illumina or Nanopore data

#### `--ref_db`

- Optional
- Type: string

Reference database in fasta file, sequence ID must be in format `SequenceName_segment#_segmentName`, for example `2021-FAV33-OS_segment1_PB2`

#### `--outdir`

- Optional
- Type: string
- Default: `./results`

The output directory where the results will be saved.

### IRMA assembly options

#### `--irma_module`

- Optional
- Type: string
- Default: `FLU-utr`

IRMA module to use for analysis.

#### `--keep_ref_deletions`

- Optional
- Type: boolean
- Default: `true`

Set "DEL_TYPE=NNN" to keep deletions to reference sequence as N characters in consensus.

### Variant Calling options

#### `--variant_caller`

- Optional
- Type: string
- Default: `clair3`

Set Variant Caller

#### `--medaka_variant_model`

- Optional
- Type: string
- Default: `r941_prom_hac_variant_g507`

Medaka model for final variant calling from phased reads

#### `--medaka_snp_model`

- Optional
- Type: string
- Default: `r941_prom_hac_snp_g507`

Medaka model for initial SNP calling from mixed reads prior to phasing

#### `--clair3_variant_model`

- Optional
- Type: string
- Default: `r941_prom_sup_g5014`

Specify built-in Clair3 variant calling model. See [Clair3 docs for options](https://github.com/HKU-BAL/Clair3#pre-trained-models)

#### `--clair3_user_variant_model`

- Optional
- Type: string
- Default: ''

Specify custom Clair3 model, e.g. [Rerio](https://github.com/nanoporetech/rerio) Clair3 model for R10 flowcells (--clair3_user_variant_model /path/to/rerio-models/r1041_e82_400bps_sup_g615_model)

#### `--minor_allele_fraction`

- Optional
- Type: number
- Default: `0.25`

Set Minor variant allele frequency/fraction

#### `--major_allele_fraction`

- Optional
- Type: number
- Default: `0.75`

Set Major variant allele frequency/fraction. Only major variant alleles are used for generating a consensus sequence

#### `--low_coverage`

- Optional
- Type: number
- Default: `10`

Low coverage depth threshold. Consensus sequence positions with less than this coverage depth will be masked with `N`

### Nanopore options

#### `--min_sample_reads`

- Optional
- Type: number
- Default: `100`

Minimum number of raw reads required per sample in order to be considered for the downstream processing steps.

### Illumina options

#### `--min_alternate_fraction`

- Optional
- Type: number
- Default: `0.5`

Minimum fraction of observations supporting an alternate allele within a single sample in the in order to evaluate the position.

### Mismatch Report options

#### `--min_aln_length`

- Optional
- Type: integer
- Default: `700`

Minimum alignment length of nucleotide BLAST results to consider for getting mismatch report from BLAST result against reference database. We normally set this number around the smallest genome segment (Segment 8)

### H/N subtyping options

Hemaglutinin and neuraminase subtype prediction options

#### `--pident_threshold`

- Optional
- Type: number
- Default: `0.85`

Minimum % identity of nucleotide BLAST results to consider for determining H/N subtypes.

#### `--max_top_blastn`

- Optional
- Type: number
- Default: `3`

Maximum of top blastn result reported

#### `--ncbi_influenza_fasta`

- Optional
- Type: string
- Default: `https://api.figshare.com/v2/file/download/41415330`

Path/URL to Zstandard compressed NCBI Influenza virus sequences FASTA file.

#### `--ncbi_influenza_metadata`

- Optional
- Type: string
- Default: `https://api.figshare.com/v2/file/download/41415333`

Path/URL to Zstandard compressed NCBI Influenza virus sequences metadata CSV file.

### Generic options

Less common options for the pipeline, typically set in a config file.

#### `--help`

- Optional
- Type: boolean

Display help text.

#### `--publish_dir_mode`

- Optional
- Type: string
- Default: `copy`

Method used to save pipeline results to output directory.

> **NOTE:** The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.

#### `--monochrome_logs`

- Optional
- Type: boolean

Do not use coloured log outputs.

> **NOTE:** Set to disable colourful command line output and live life in monochrome.

#### `--tracedir`

- Optional
- Type: string
- Default: `${params.outdir}/pipeline_info`

Directory to keep pipeline Nextflow logs and reports.

#### `--singularity_pull_docker_container`

- Optional
- Type: boolean

Instead of directly downloading Singularity images for use with Singularity, force the workflow to pull and convert Docker containers instead.

#### `--validate_params`

- Optional
- Type: boolean
- Default: `True`

Boolean whether to validate parameters against the schema at runtime

#### `--show_hidden_params`

- Optional
- Type: boolean

Show all params when using `--help`

> **NOTE:** By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.

### Slurm scheduler options

#### `--slurm_queue`

- Optional
- Type: string

Slurm queue/partition to submit pipeline jobs to

#### `--slurm_queue_size`

- Optional
- Type: integer
- Default: `100`

Slurm queue size. Max number of jobs to queue at once.

### Max job request options

Set the top limit for requested resources for any single job.

#### `--max_cpus`

- Optional
- Type: integer
- Default: `16`

Maximum number of CPUs that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`

#### `--max_memory`

- Optional
- Type: string
- Default: `128.GB`

Maximum amount of memory that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`

#### `--max_time`

- Optional
- Type: string
- Default: `240.h`

Maximum amount of time that can be requested for any single job.

> **NOTE:** Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test_illumina` and `test_nanopore`
  - Profiles for testing the Illumina and Nanopore workflows, respectively, with complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
