# CFIA-NCFAD/nf-flu - Influenza A and B Virus Genome Assembly Nextflow Workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15093332.svg)](https://doi.org/10.5281/zenodo.15093332)
[![CI](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml/badge.svg)](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with podman](https://img.shields.io/badge/run%20with-podman-1d355c.svg?labelColor=000000)](https://podman.io/)

## Introduction

**nf-flu** is a [Nextflow][] bioinformatics analysis pipeline for assembly and analysis of Influenza A and B viruses from Illumina or Nanopore sequencing data or previously assembled FASTA sequences.
Since Influenza has a segmented genome consisting of 8 gene segments, the pipeline will automatically select the top matching reference sequence from NCBI for each gene segment based on [IRMA][] assembly and nucleotide [BLAST][] against all Influenza sequences from NCBI.
Users can also provide their own reference sequences to include in the top reference sequence selection process.
After reference sequence selection, the pipeline performs read mapping to each reference sequence, variant calling and depth-masked consensus sequence generation.

> **Note:** The officially supported version of the pipeline is [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu). If you have issues with using the pipeline, please [create an issue](https://github.com/CFIA-NCFAD/nf-flu/issues/new/choose) on [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu) repo.

## Pipeline summary

1. Download latest NCBI Influenza virus sequences and metadata (see [docs](docs/update_seqs_db.md) for more details).
2. Merge reads of re-sequenced samples ([`cat`](http://www.linfo.org/cat.html)) (if needed).
3. Assembly of Influenza gene segments with [IRMA][] using the built-in FLU module
4. Nucleotide [BLAST][] search against [NCBI Influenza DB][] sequences
5. H/N subtype prediction and Excel XLSX report generation based on BLAST results.
6. Automatically select top match reference sequences for segments
7. Read mapping, variant calling and consensus sequence generation for each segment against top reference sequence based on BLAST results.
8. Annotation of consensus sequences with [VADR][]
9. [FluMut][] detection of molecular markers and mutation in Influenza A(H5N1) viruses.
10. [GenoFLU][] genotyping of North American H5 viruses. [Genin2][] genotyping using information from clade 2.3.4.4b H5Nx viruses collected in Europe since October 2020.
11. HA cleavage site prediction and classification
12. [Nextclade][] clade assignment, mutation calling and sequence quality checks.
13. [MultiQC][] report generation.

![nf-flu workflow](assets/nf-flu-pipeline-diagram.svg)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`; latest stable release recommended!).
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Apptainer`][], [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort)_
3. Download the pipeline and test it on a minimal dataset with a single command:

    For Illumina workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_illumina,<docker/apptainer/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    For Nanopore workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_nanopore,<docker/apptainer/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    > * If you are using `apptainer`/`singularity` then the pipeline will auto-detect this and attempt to download the Apptainer/Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Apptainer/Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Run your own analysis

    * [Optional] Generate an input samplesheet from a directory containing Illumina FASTQ files (e.g. `/path/to/illumina_run/Data/Intensities/Basecalls/`) with the included Python script [`fastq_dir_to_samplesheet.py`](bin/fastq_dir_to_samplesheet.py) **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```bash
        python ~/.nextflow/assets/CFIA-NCFAD/nf-flu/bin/fastq_dir_to_samplesheet.py \
          -i /path/to/illumina_run/Data/Intensities/Basecalls/ \
          -o samplesheet.csv
        ```

    * Typical command for Illumina sequencing data

        ```bash
        nextflow run CFIA-NCFAD/nf-flu \
          --input samplesheet.csv \
          --platform illumina \
          --profile <docker/apptainer/singularity/podman/shifter/charliecloud/conda>
        ```

    * Typical command for Nanopore sequencing data

      ```bash
      nextflow run CFIA-NCFAD/nf-flu \
        --input samplesheet.csv \
        --platform nanopore \
        --profile <docker/apptainer/singularity/conda>
      ```

    * Run analysis on FASTA files within a directory

      ```bash
      nextflow run CFIA-NCFAD/nf-flu \
        --input /path/to/fasta_files/ \
        --platform assemblies \
        --profile <docker/apptainer/singularity/conda>
      ```

## Documentation

The nf-flu pipeline comes with:

* [Usage](docs/usage.md) and
* [Output](docs/output.md) documentation.

## Resources and References

### [BcfTools][] and [Samtools][]

```text
Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H., 2021. Twelve years of SAMtools and BCFtools. Gigascience 10, giab008.
https://doi.org/10.1093/gigascience/giab008
```

### [BLAST][] Basic Local Alignment Search Tool

```text
Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J., 1990. Basic local alignment search tool. J. Mol. Biol. 215, 403–410.
https://doi.org/10.1016/S0022-2836(05)80360-2
```

```text
Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L., 2009. BLAST+: architecture and applications. BMC Bioinformatics 10, 421.
https://doi.org/10.1186/1471-2105-10-421
```

### [Clair3][]

```text
Zheng, Z., Li, S., Su, J., Leung, A.W.-S., Lam, T.-W., Luo, R., 2022. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. Nat Comput Sci 2, 797–803.
https://doi.org/10.1038/s43588-022-00387-x
```

### [FluMut][]

[FluMut][] is used to "search for molecular markers with potential impact on the biological characteristics of Influenza A viruses of the A(H5N1) subtype, starting from complete or partial nucleotide genome sequences".

```text
Edoardo Giussani, Alessandro Sartori, Angela Salomoni, Lara Cavicchio, Cristian de Battisti, Ambra Pastori, Maria Varotto, Bianca Zecchin, Joseph Hughes, Isabella Monne, Alice Fusaro. FluMut: a tool for mutation surveillance in highly pathogenic H5N1 genomes. Virus Evolution, Volume 11, Issue 1, 2025, veaf011.
https://doi.org/10.1093/ve/veaf011
```

### [Freebayes][]

[Freebayes][] is used for variant calling.

```text
Garrison, E., Marth, G., 2012. Haplotype-based variant detection from short-read sequencing. arXiv:1207.3907 [q-bio].
https://doi.org/10.48550/arXiv.1207.3907
```

### [Genin2][]

[Genin2][] "predicts genotypes for clade 2.3.4.4b H5Nx viruses collected in Europe since October 2020".

```text
https://github.com/izsvenezie-virology/genin2
```

> **NOTE:** Authors [suggest citing the GitHub repo](https://github.com/izsvenezie-virology/genin2?tab=readme-ov-file#cite-genin2) at this time.

### [GenoFLU][]

[GenoFLU][] "identifies the genotype of North American H5 2.3.4.4b viruses as well as providing information on individual segments when a sequence does not belong to a defined genotype".

```text
Youk S, Torchetti MK, Lantz K, Lenoch JB, Killian ML, Leyson C, Bevins SN, Dilione K, Ip HS, Stallknecht DE, Poulson RL, Suarez DL, Swayne DE, Pantin-Jackwood MJ. 
H5N1 highly pathogenic avian influenza clade 2.3.4.4b in wild and domestic birds: Introductions into the United States and reassortments, December 2021-April 2022. Virology. 2023 Oct;587:109860. doi: 10.1016/j.virol.2023.109860. Epub 2023 Aug 2. PMID: 37572517.
```

> **NOTE:** The authors recommend citing the [Youk et al paper](https://doi.org/10.1016/j.virol.2023.109860) according to [GenoFLU issue #10](https://github.com/USDA-VS/GenoFLU/issues/10).

### [IRMA][] Iterative Refinement Meta-Assembler

```text
Shepard, S.S., Meno, S., Bahl, J., Wilson, M.M., Barnes, J., Neuhaus, E., 2016. Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler. BMC Genomics 17, 708.
https://doi.org/10.1186/s12864-016-3030-6
```

### [Medaka][]

[Medaka][] is deprecated in favour of [Clair3][] for variant calling of Nanopore data.

### [Minimap2][]

[Minimap2][] is used for rapid and accurate read alignment to reference sequences.

```text
Li, H., 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34, 3094–3100.
https://doi.org/10.1093/bioinformatics/bty191
```

### [Mosdepth][]

[Mosdepth][] is used for rapid sequencing coverage calculation and summary statistics.

```text
Pedersen, B.S., Quinlan, A.R., 2017. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics 34, 867–868.
https://doi.org/10.1093/bioinformatics/btx699
```

### [MultiQC][]

[MultiQC][] is used for generation of a single report for multiple tools.

```text
Ewels, P., Magnusson, M., Lundin, S., Käller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047–3048.
https://doi.org/10.1093/bioinformatics/btw354
```

### [NCBI Influenza Virus Resource][]

**nf-flu** relies on publicly available Influenza sequence data from NCBI available at the [NCBI Influenza Virus Resource][], which is downloaded from the [FTP site](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/).

NCBI Influenza Virus Resource:

```text
Bao, Y., Bolotov, P., Dernovoy, D., Kiryutin, B., Zaslavsky, L., Tatusova, T., Ostell, J., Lipman, D., 2008. The influenza virus resource at the National Center for Biotechnology Information. J Virol 82, 596–601.
https://doi.org/10.1128/JVI.02005-07
```

NCBI Influenza Virus Sequence Annotation Tool:

```text
Bao, Y., Bolotov, P., Dernovoy, D., Kiryutin, B., Tatusova, T., 2007. FLAN: a web server for influenza virus genome annotation. Nucleic Acids Res 35, W280-284.
https://doi.org/10.1093/nar/gkm354
```

### [Nextclade][]

**nf-flu** performs Nextclade clade assignment, mutation calling and sequence quality checks of assembled Influenza sequences against 30 Nextclade datasets for different subtypes and lineages of Influenza A and B virus.

```text
Aksamentov, I., Roemer, C., Hodcroft, E. B., & Neher, R. A., (2021). Nextclade: clade assignment, mutation calling and quality control for viral genomes. Journal of Open Source Software, 6(67), 3773,
https://doi.org/10.21105/joss.03773
```

### [Nextflow][]

**nf-flu** is implemented in [Nextflow][].

```text
Tommaso, P.D., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., Notredame, C., 2017. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319.
https://doi.org/10.1038/nbt.3820
```

### [nf-core][]

[nf-core][] is a great resource for building robust and reproducible bioinformatics pipelines.

```text
Ewels, P.A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M.U., Di Tommaso, P., Nahnsen, S., 2020. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38, 276–278.
https://doi.org/10.1038/s41587-020-0439-x
```

### [seqtk][]

[seqtk][] is used for rapid manipulation of FASTA/Q files. Available from GitHub at [lh3/seqtk](https://github.com/lh3/seqtk)

### [VADR][]

[VADR][] is used for annotation of Influenza virus sequences.

```text
Alejandro A Schäffer, Eneida L Hatcher, Linda Yankie, Lara Shonkwiler, J Rodney Brister, Ilene Karsch-Mizrachi, Eric P Nawrocki; VADR: validation and annotation of virus sequence submissions to GenBank. BMC Bioinformatics 21, 211 (2020).
https://doi.org/10.1186/s12859-020-3537-3
```

### [table2asn][]

[table2asn][] is used for converting the [VADR][] Feature Table format output to Genbank format to help with conversion to other formats such as FASTA and GFF.

## Contributors

* [Peter Kruczkiewicz](https://github.com/peterk87) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - lead developer
* [Hai Nguyen](https://github.com/nhhaidee) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - Nanopore workflow
* [Abdallah Meknas](https://github.com/ameknas-phac) (Influenza, Respiratory Viruses, and Coronavirus Section (IRVC), Public Health Agency of Canada (PHAC)) - expansion of the Illumina workflow
* [Cass Erdelyan](https://github.com/cerdelyan/) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - development, testing and valuable feedback

## Credits

* [nf-core](https://nf-co.re) project for establishing Nextflow workflow development best-practices, [nf-core tools](https://nf-co.re/tools-docs/) and [nf-core modules](https://github.com/nf-core/modules)
* [nf-core/viralrecon](https://github.com/nf-core/viralrecon) for inspiration and setting a high standard for viral sequence data analysis pipelines
* [Conda](https://docs.conda.io/projects/conda/en/latest/) and [Bioconda](https://bioconda.github.io/) project for making it easy to install, distribute and use bioinformatics software.
* [Biocontainers](https://biocontainers.pro/) for automatic creation of [Docker] and [Apptainer]/[Singularity] containers for bioinformatics software in [Bioconda]

[Apptainer]: https://apptainer.org/
[BcfTools]: https://samtools.github.io/bcftools/
[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[Clair3]: https://github.com/HKU-BAL/Clair3
[Docker]: https://www.docker.com/
[FluMut]: https://github.com/izsvenezie-virology/FluMut
[Freebayes]: https://github.com/freebayes/freebayes
[Genin2]: https://github.com/izsvenezie-virology/genin2
[GenoFLU]: https://github.com/USDA-VS/GenoFLU
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[Medaka]: https://github.com/nanoporetech/medaka
[Minimap2]: https://github.com/lh3/minimap2/
[Mosdepth]: https://github.com/brentp/mosdepth
[MultiQC]: https://multiqc.info/
[NCBI Influenza DB]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[NCBI Influenza Virus Resource]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[Nextclade]: https://clades.nextstrain.org/
[Nextflow]: https://www.nextflow.io/
[nf-core]: https://nf-co.re/
[Samtools]: https://www.htslib.org/
[seqtk]: https://github.com/lh3/seqtk
[Singularity]: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
[table2asn]: https://www.ncbi.nlm.nih.gov/genbank/table2asn/
[VADR]: https://github.com/ncbi/vadr
