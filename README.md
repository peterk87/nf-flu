# CFIA-NCFAD/nf-flu - Influenza A and B Virus Genome Assembly Nextflow Workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13892044.svg)](https://doi.org/10.5281/zenodo.13892044)
[![CI](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml/badge.svg)](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-flu** is a [Nextflow][] bioinformatics analysis pipeline for assembly and H/N subtyping of Influenza A and B viruses from Illumina or Nanopore sequencing data.
Since Influenza has a segmented genome consisting of 8 gene segments, the pipeline will automatically select the top matching reference sequence from NCBI for each gene segment based on [IRMA][] assembly and nucleotide [BLAST][] against all Influenza sequences from NCBI.
Users can also provide their own reference sequences to include in the top reference sequence selection process.
After reference sequence selection, the pipeline performs read mapping to each reference sequence, variant calling and depth-masked consensus sequence generation.

## Pipeline summary

1. Download latest [NCBI Orthomyxoviridae sequences](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=11308&lvl=3&keep=1&srchmode=1&unlock) and metadata (parsed from [NCBI Viruses FTP data](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/)).
2. Merge reads of re-sequenced samples ([`cat`](http://www.linfo.org/cat.html)) (if needed).
3. Assembly of Influenza gene segments with [IRMA][] using the built-in FLU module
4. Nucleotide [BLAST][] search against [NCBI Influenza DB][] sequences
5. H/N subtype prediction and Excel XLSX report generation based on BLAST results.
6. Automatically select top match reference sequences for segments
7. Read mapping, variant calling and consensus sequence generation for each segment against top reference sequence based on BLAST results.
8. Annotation of consensus sequences with [VADR][]
9. [MultiQC][] report generation.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`).
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort)_
3. Download the pipeline and test it on a minimal dataset with a single command:

    For Illumina workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_illumina,<docker/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    For Nanopore workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_nanopore,<docker/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    > * If you are using `singularity` then the pipeline will auto-detect this and attempt to download the Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Run your own analysis

    * [Optional] Generate an input samplesheet from a directory containing Illumina FASTQ files (e.g. `/path/to/illumina_run/Data/Intensities/Basecalls/`) with the included Python script [`fastq_dir_to_samplesheet.py`](bin/fastq_dir_to_samplesheet.py) **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```bash
        python ~/.nextflow/assets/CFIA-NCFAD/nf-flu/bin/fastq_dir_to_samplesheet.py \
          -i /path/to/illumina_run/Data/Intensities/Basecalls/ \
          -o samplesheet.csv
        ```

    * Typical command for Illumina Platform

        ```bash
        nextflow run CFIA-NCFAD/nf-flu \
          --input samplesheet.csv \
          --platform illumina \
          --profile <docker/singularity/podman/shifter/charliecloud/conda>
        ```

    * Typical command for Nanopore Platform

      ```bash
      nextflow run CFIA-NCFAD/nf-flu \
        --input samplesheet.csv \
        --platform nanopore \
        --profile <docker/singularity/conda>
      ```

## Documentation

The nf-flu pipeline comes with:

* [Usage](docs/usage.md) and
* [Output](docs/output.md) documentation.

## Resources and References

### [BcfTools][] and [Samtools][]

```text
Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H., 2021. Twelve years of SAMtools and BCFtools. Gigascience 10, giab008. https://doi.org/10.1093/gigascience/giab008
```

### [BLAST][] Basic Local Alignment Search Tool

```text
Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J., 1990. Basic local alignment search tool. J. Mol. Biol. 215, 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2
```

```text
Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., Madden, T.L., 2009. BLAST+: architecture and applications. BMC Bioinformatics 10, 421. https://doi.org/10.1186/1471-2105-10-421
```

### [Clair3][]

```text
Zheng, Z., Li, S., Su, J., Leung, A.W.-S., Lam, T.-W., Luo, R., 2022. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. Nat Comput Sci 2, 797–803. https://doi.org/10.1038/s43588-022-00387-x
```

### [IRMA][] Iterative Refinement Meta-Assembler

```text
Shepard, S.S., Meno, S., Bahl, J., Wilson, M.M., Barnes, J., Neuhaus, E., 2016. Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler. BMC Genomics 17, 708. https://doi.org/10.1186/s12864-016-3030-6
```

### [Medaka][]

[Medaka][] is deprecated in favour of [Clair3][] for variant calling of Nanopore data.

### [Minimap2][]

[Minimap2][] is used for rapid and accurate read alignment to reference sequences.

```text
Li, H., 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34, 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
```

### [Mosdepth][]

[Mosdepth][] is used for rapid sequencing coverage calculation and summary statistics.

```text
Pedersen, B.S., Quinlan, A.R., 2017. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics 34, 867–868. https://doi.org/10.1093/bioinformatics/btx699
```

### [MultiQC][]

[MultiQC][] is used for generation of a single report for multiple tools.

```text
Ewels, P., Magnusson, M., Lundin, S., Käller, M., 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
```

### [NCBI Influenza Virus Resource][]

**nf-flu** relies on publicly available Influenza sequence data from NCBI available at the [NCBI Influenza Virus Resource][], which is downloaded from the [FTP site](https://ftp.ncbi.nih.gov/genomes/INFLUENZA/).

NCBI Influenza Virus Resource:

```text
Bao, Y., Bolotov, P., Dernovoy, D., Kiryutin, B., Zaslavsky, L., Tatusova, T., Ostell, J., Lipman, D., 2008. The influenza virus resource at the National Center for Biotechnology Information. J Virol 82, 596–601. https://doi.org/10.1128/JVI.02005-07
```

NCBI Influenza Virus Sequence Annotation Tool:

```text
Bao, Y., Bolotov, P., Dernovoy, D., Kiryutin, B., Tatusova, T., 2007. FLAN: a web server for influenza virus genome annotation. Nucleic Acids Res 35, W280-284. https://doi.org/10.1093/nar/gkm354
```

### [Nextflow][]

**nf-flu** is implemented in [Nextflow][].

```text
Tommaso, P.D., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., Notredame, C., 2017. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319. https://doi.org/10.1038/nbt.3820
```

### [nf-core][]

[nf-core][] is a great resource for building robust and reproducible bioinformatics pipelines.

```text
Ewels, P.A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M.U., Di Tommaso, P., Nahnsen, S., 2020. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38, 276–278. https://doi.org/10.1038/s41587-020-0439-x
```

### [seqtk][]

[seqtk][] is used for rapid manipulation of FASTA/Q files. Available from GitHub at [lh3/seqtk](https://github.com/lh3/seqtk)

### [VADR][]

[VADR][] is used for annotation of Influenza virus sequences.

```text
Alejandro A Schäffer, Eneida L Hatcher, Linda Yankie, Lara Shonkwiler, J Rodney Brister, Ilene Karsch-Mizrachi, Eric P Nawrocki; VADR: validation and annotation of virus sequence submissions to GenBank. BMC Bioinformatics 21, 211 (2020). https://doi.org/10.1186/s12859-020-3537-3
```

### [table2asn][]

[table2asn][] is used for converting the [VADR][] Feature Table format output to Genbank format to help with conversion to other formats such as FASTA and GFF.

## Contributors

* [Peter Kruczkiewicz](https://github.com/peterk87) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - lead developer
* [Hai Nguyen](https://github.com/nhhaidee) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - Nanopore workflow
* [Abdallah Meknas](https://github.com/ameknas-phac) (Influenza, Respiratory Viruses, and Coronavirus Section (IRVC), Public Health Agency of Canada (PHAC)) - expansion of the Illumina workflow

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
[NCBI Influenza Virus Resource]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[BcfTools]: https://samtools.github.io/bcftools/
[Samtools]: https://www.htslib.org/
[nf-core]: https://nf-co.re/
[Minimap2]: https://github.com/lh3/minimap2/
[Clair3]: https://github.com/HKU-BAL/Clair3
[Medaka]: https://github.com/nanoporetech/medaka
[Mosdepth]: https://github.com/brentp/mosdepth
[seqtk]: https://github.com/lh3/seqtk
[MultiQC]: https://multiqc.info/
[VADR]: https://github.com/ncbi/vadr
[table2asn]: https://www.ncbi.nlm.nih.gov/genbank/table2asn/
