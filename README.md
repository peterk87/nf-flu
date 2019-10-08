# Influenza Genome Analysis Nextflow Workflow

## Introduction

This repo contains a [Nextflow] workflow for the [IRMA] assembly and H/N subtyping by nucleotide [BLAST] against the [NCBI Influenza DB] of Illumina sequenced influenza virus genomes.

## Pipeline Steps

| Step                                                | Main program/s                      |
|-----------------------------------------------------|-------------------------------------|
| IRMA iterative assembly of genome segments          | IRMA                                |
| H/N subtyping                                       | BLAST+                              |


## Getting Started

1. Install [Nextflow][] (Java 8 or later must be installed; Conda provides an easy way to install Nextflow with all dependencies required)

```bash
curl -s https://get.nextflow.io | bash
# install nextflow globally for all users
sudo mv nextflow /usr/local/bin/
# or locally to your user
mkdir -p ~/bin
mv nextflow ~/bin/
```

2. Install [Singularity][] 

3. Run this workflow

```bash
nextflow run peterk87/nf-iav-illumina --help
```

## Usage

### Run on directory of Illumina paired-end reads

```
nextflow run peterk87/nf-iav-illumina --reads_dir /path/to/reads --outdir /path/to/results/outdir
```

### Run on a cluster with SLURM

```
nextflow run peterk87/nf-iav-illumina --reads_dir /path/to/reads -profile slurm --slurm_queue SLURM_QUEUE_NAME --outdir /path/to/results/outdir
```


### Show help info

```
nextflow run peterk87/nf-iav-illumina --help
```

```
===========================================
peterk87/nf-iav-illumina  ~  version 1.1.0
===========================================

Git info: null - null [null]

Usage:
The typical command for running the pipeline is as follows:
nextflow run peterk87/nf-iav-illumina --reads_dir "/path/to/reads/*R{1,2}*.fastq.gz" --outdir /path/to/results
Mandatory arguments:
  --reads                   Illumina FASTQ reads
Other options:
  --outdir                  The output directory where the results will be saved
  -w/--work-dir             The temporary directory where intermediate data will be saved
  --slurm_queue             Name of SLURM queue to run workflow on (default "")
  -profile                  Configuration profile to use. [standard, slurm] (default 'standard')
```


## Results

This workflow outputs results in an output directory with the following structure:

```
<results output directory>
├── consensus
│   ├── Sample1
│   │   ├── Sample1_1.fa
│   │   ├── ...
│   │   └── Sample1_8.fa
│   ├── Sample2
│   │   ├── Sample2_1.fa
│   │   ├── ...
│   │   └── Sample2_8.fa
├── irma
│   ├── Sample1
│   │   └── Sample1-read-counts.tsv
│   ├── Sample2
│   │   └── Sample2-read-counts.tsv
│   ├── Sample1.tar.gz
│   ├── Sample2.tar.gz
├── irma-consensus-report.tsv
├── pipeline_info
│   ├── nf-iav-illumina_DAG.svg
│   ├── nf-iav-illumina_report.html
│   ├── nf-iav-illumina_timeline.html
│   └── nf-iav-illumina_trace.txt
└── subtyping_report.xlsx
```

- `consensus` contains [IRMA] generated consensus sequences for each sample for each of the genome segments (if available)
- `irma` contains a `.tar.gz` of the [IRMA] output directory for each sample
- `irma-consensus-report.tsv` is a tab-delimited table of the number of [IRMA] consensus genome segments obtained for each sample
- `subtyping_report.xlsx` is summary report of H/N subtyping by [BLAST] containing subtyping results for all input samples



### Workflow Execution Graph

![Workflow execution directed-acyclic graph](dag.svg)




[NCBI Influenza DB]: ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/
[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[Nextflow]: https://www.nextflow.io/
[Singularity]: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
