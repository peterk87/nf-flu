# Influenza Genome Analysis Nextflow Workflow

[![Build Status](https://dev.azure.com/peterkruczkiewicz/nf-iav-illumina/_apis/build/status/peterk87.nf-iav-illumina?branchName=master)](https://dev.azure.com/peterkruczkiewicz/nf-iav-illumina/_build/latest?definitionId=1&branchName=master)

## Introduction

This repo contains a [Nextflow][] workflow for the [IRMA][] assembly and H/N subtyping by nucleotide [BLAST][] against the [NCBI Influenza DB][].

## Pipeline Steps

| Step                                                | Main program                        |
|-----------------------------------------------------|-------------------------------------|
| IRMA iterative assembly of genome segments          | IRMA                                |
| H/N subtyping                                       | BLAST+                              |

Singularity containers used in this workflow:

| Singularity container | Hosted At | Description |
|-----------------------|-----------|-------------|
| [peterk87/nf-iav-illumina](https://singularity-hub.org/collections/3633) | [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3633) | Singularity container for [BLAST][]+ and [BLAST][] result parsing processes |
| [peterk87/irma-singularity](https://singularity-hub.org/collections/2385) | [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2385) | Singularity container for [IRMA][] |

## Getting Started

1. Install [Nextflow][] (Java 8 or later must be installed; Conda provides an easy way to install Nextflow with all dependencies required)

Conda install:

```bash
# Setup Conda channels for Bioconda and Conda Forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Install Nextflow with Conda
conda install nextflow
```

Manual install:

```bash
curl -s https://get.nextflow.io | bash
# install nextflow globally for all users
sudo mv nextflow /usr/local/bin/
# or locally to your user
mkdir -p ~/bin
mv nextflow ~/bin/
```

2. Install [Singularity][] (required for running workflow both locally and in a high-performance computing cluster)

3. Try to show help info for this workflow

```bash
nextflow run peterk87/nf-iav-illumina --help
```

4. Run on some test data

Download reads for Influenza A virus (A/England/195/2009(H1N1)) (ERR3338653; https://www.ncbi.nlm.nih.gov/sra/ERX3363362[accn]) with [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) and this Nextflow workflow:

```bash
mkdir -p test/reads && cd test/reads
fasterq-dump ERR3338653
# change FASTQ headers to make them compatible with IRMA/LABEL, otherwise, 
# empty consensus sequences are produced and the only log message produced is
# "Irregular header for fastQ read pairs"
mkdir original && mv ERR3338653*.fastq original/
cat original/ERR3338653_1.fastq | perl -pe 's/^@([\w\.]+) (\S+) (\S+)/@\2 1:N:0:51/' > ERR3338653_1.fastq
cat original/ERR3338653_2.fastq | perl -pe 's/^@([\w\.]+) (\S+) (\S+)/@\2 2:N:0:51/' > ERR3338653_2.fastq
cd ..
# Run the workflow!
nextflow run peterk87/nf-iav-illumina --reads "reads/*_{1,2}.fastq" --outdir results
```

You should find the subtype to be H1N1 for ERR3338653.

## Usage

It's a good idea to run the latest version of the workflow if possible. 
You can update the workflow with:

```bash
nextflow pull peterk87/nf-iav-illumina
```

This will pull the latest version of the workflow from GitHub.

### Run locally on directory of Illumina paired-end reads

*NOTE: Please ensure that your workstation has [Singularity][] installed! Run `which singularity` and `singularity --version` to ensure that a recent version of Singularity has been installed.*

```bash
nextflow run peterk87/nf-iav-illumina \
  --reads "reads/*_R{1,2}*.fastq.gz" \
  --outdir /path/to/results/outdir
```

### Run on a high-performance computing (HPC) cluster with SLURM

*NOTE: Please ensure that your cluster has [Singularity][] installed! Run `which singularity` and `singularity --version` to ensure that a recent version of Singularity has been installed.*

```bash
nextflow run peterk87/nf-iav-illumina \
  --reads "reads/*_R{1,2}*.fastq.gz" \
  --outdir /path/to/results/outdir \
  -profile slurm \
  --slurm_queue SLURM_QUEUE_NAME 
```

Replace `SLURM_QUEUE_NAME` with the name of the queue or partition you can submit to on your cluster.
Contact your cluster tech support to find out which queue you can submit to on the cluster. 

Show SLURM info about your cluster with the `sinfo` command and `squeue` to show how busy the cluster is. 

### Show help info

```bash
nextflow run peterk87/nf-iav-illumina --help
```

```
===========================================
peterk87/nf-iav-illumina  ~  version 1.1.0
===========================================

  Git info: null - null [null]

Usage:

The typical command for running the pipeline is as follows:

  nextflow run peterk87/nf-iav-illumina --reads "reads/*R{1,2}*.fastq.gz" --outdir results

NOTE: Please ensure you have Singularity installed prior to running this workflow. (https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)

Mandatory Options:
  --reads           Input paired-end Illumina FASTQ reads; quotes required! (default: "reads/*R{1,2}*.fastq.gz")
Other Options:
  --outdir          Output results directory (default: results)
  -w/--work-dir     The temporary directory where intermediate data will be saved (default: work)
  -profile          Configuration profile to use [standard, slurm] (default: "standard")
Cluster Options:
  --slurm_queue     Name of SLURM queue to run workflow on; use with -profile slurm
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

### Subtyping report spreadsheet

`subtyping_report.xlsx` is summary report of H/N subtyping by [BLAST] containing subtyping results for all input samples. 

The subtyping report spreadsheet contains four sheets:

- **1_Subtype Predictions**: H/N subtype prediction results for each sample along with top matching [NCBI Influenza DB][] segment for the H and N segments
- **2_Top Segment Matches**: Top 3 [NCBI Influenza DB][] sequence matches for each segment of each sample along with BLASTN hit values and reference sequence metadata. 
- **3_H Segment Results**: Top H subtype prediction, BLASTN results and top matching sequence metadata for each sample.
- **4_N Segment Results**: Top N subtype prediction, BLASTN results and top matching sequence metadata for each sample.

#### Sheet: 1_Subtype Predictions

This sheet contains the H/N subtype prediction results for each sample along with top matching [NCBI Influenza DB][] segment for the H and N segments

| Field | Description | Example |
|-------|-------------|---------|
| Sample | Sample name | ERR3338653 | 
| Subtype Prediction | H/N subtype prediction based on BLAST+ against the [NCBI Influenza DB][]. If a type could not be assigned to either H or N segment or both, then the subtype prediction will be missing the H or N value or if both the H and N cannot be assigned then the subtype prediction will be null or an empty cell value  | H1N1 |
| H: top match accession | NCBI accession of top matching [NCBI Influenza DB] sequence for the H segment | CY147779 |
| H: type prediction | H subtype prediction number. Value is a number. | 1 |
| H: top match virus name | Top matching sequence virus name | Influenza A virus (A/Mexico/24036/2009(H1N1)) |
| H: NCBI Influenza DB subtype match proportion | Proportion of BLAST matches that support the H subtype prediction. This value is a decimal number where 1.0 indicates 100% of matches support the subtype prediction. | 0.9980057896 | 
| N: top match accession | NCBI accession of top matching [NCBI Influenza DB] sequence for the N segment | MN371610 |
| N: type prediction | N subtype prediction. Value is a number. | 1 |
| N: top match virus name  | Top matching sequence virus name | Influenza A virus (A/California/04/2009) |
| N: NCBI Influenza DB subtype match proportion | Proportion of BLAST matches that support the N subtype prediction. This value is a decimal number where 1.0 indicates 100% of matches support the subtype prediction. | 0.9993240503 |



#### Sheet: 2_Top Segment Matches

This sheet contains the top 3 [NCBI Influenza DB][] sequence matches for each segment of each sample along with BLASTN hit values and reference sequence metadata. 



| Field | Description | Example |
|-------|-------------|---------|
| Sample | Sample name | ERR3338653 |
| Sample Genome Segment Number | Influenza genome segment number | 4 | 
| Reference NCBI Accession | Matching sequence NCBI accession | CY147779 |
| Reference Subtype | Matching sequence subtype | H1N1 |
| BLASTN Percent Identity | BLASTN percent identity | 99.941 |
| BLASTN Alignment Length | BLASTN alignment length | 1701 |
| BLASTN Mismatches | BLASTN number of mismatches | 1 |
| BLASTN Gaps | BLASTN number of gaps | 0 |
| BLASTN Sample Start Index | Sample sequence alignment start index | 1 |
| BLASTN Sample End Index | Sample sequence alignment end index | 1701 |
| BLASTN Reference Start Index | Matching reference sequence start index | 1 |
| BLASTN Reference End Index | Matching reference sequence end index | 1701 |
| BLASTN E-value | BLASTN alignment e-value | 0 |
| BLASTN Bitscore | BLASTN alignment bitscore | 3136 |
| Sample Sequence Length | Length of sample sequence segment | 1701 |
| Reference Sequence Length | Length of matching reference sequence segment | 1701 |
| Sample Sequence Coverage of Reference Sequence | Sample segment sequence coverage of reference sequence  | 100 |
| Reference Sequence ID | Matching reference sequence identifier | gi|525340040|gb|CY147779|Influenza A virus (A/Mexico/24036/2009(H1N1)) hemagglutinin (HA) gene, complete cds |
| Reference Genome Segment Number | Reference sequence segment number | 4 |
| Reference Virus Name | Reference sequence virus name | Influenza A virus (A/Mexico/24036/2009(H1N1)) |
| Reference Host | Reference sequence host organism | Human |
| Reference Country | Reference sequence country of isolation | Mexico |
| Reference Collection Date | Reference sequence date of collection | 2009/04/27 |
| Reference Patient Age | Reference sequence patient age | 152 |
| Reference Patient Gender | Reference sequence patient gender | Female |
| Reference Group ID | [NCBI Influenza DB][] reference sequence internal genome ID | 1018714 |


#### Sheets: "3_H Segment Results" and "4_N Segment Results"

These sheets ("3_H Segment Results" and "4_N Segment Results") contain the subtype prediction, BLASTN results and top matching sequence metadata for each sample.  

Below are shown the fields for the "3_H Segment Results" sheet. The fields are nearly identical for the "4_N Segment Results" sheet except "N:" instead of "H:" in the field names.

| Field | Description | Example |
|-------|-------------|---------|
| Sample | Sample name | ERR3338653 |
| Subtype Prediction | Overall subtype prediction | H1N1 |
| H: NCBI Influenza DB subtype match proportion | Proportion of BLAST matches that support the N subtype prediction. This value is a decimal number where 1.0 indicates 100% of matches support the subtype prediction. | 0.9980057896 |
| H: NCBI Influenza DB subtype match count | Number of reference sequences that have a  BLASTN match to the sequence of this sample and have the same H subtype as the top prediction. | 31028 |
| H: NCBI Influenza DB total count | Total number of reference sequences that have a  BLASTN match to the sequence of this sample. | 31090 |
| H: top match BLASTN % identity | BLASTN alignment percent identity | 99.941 |
| H: top match BLASTN alignment length | BLASTN alignment length | 1701 |
| H: top match BLASTN mismatches | BLASTN alignment mismatches | 1 |
| H: top match BLASTN gaps | BLASTN alignment gaps | 0 |
| H: top match BLASTN bitscore | BLASTN alignment bitscore | 3136 |
| H: sample segment length | Sample sequence length | 1701 |
| H: top match sequence length | Reference sequence length | 1701 |
| H: top match accession | Matching sequence NCBI accession | CY147779 |
| H: top match virus name | Reference sequence virus name | Influenza A virus (A/Mexico/24036/2009(H1N1)) |
| H: top match host | Reference sequence host organism | Human |
| H: top match country | Reference sequence country of isolation | Mexico |
| H: top match collection date | Reference sequence date of collection | 2009/04/27 |
| H: type prediction | H segment subtype prediction | 1 |


## Workflow Execution Graph

![Workflow execution directed-acyclic graph](dag.svg)

## Resources

- NCBI Influenza FTP site: ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/
- [IRMA][] Iterative Refinement Meta-Assembler
  - [IRMA Publication](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6)



[NCBI Influenza DB]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[Nextflow]: https://www.nextflow.io/
[Singularity]: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
