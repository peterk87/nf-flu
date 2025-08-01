# nf-flu: Output

The output produced by the CFIA-NCFAD/nf-flu pipeline is described here.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

- [IRMA](#irma)
- [BLAST analysis](#blast-analysis)
- [Coverage Plots](#coverage-plots)
- [Assembled Consensus Sequences](#assembled-consensus-sequences)
- [Mismatch Report](#mismatch-report)
- [Reference Sequences](#reference-sequences)
- [Variant Calling](#variant-calling)
- [H/N Subtyping](#hn-subtyping)
- [Annotation](#annotation)
- [FluMut](#flumut)
- [Genin2](#genin2)
- [GenoFLU](#genoflu)
- [Cleavage Site Prediction](#cleavage-site-prediction)
- [Nextclade](#nextclade)
- [Pipeline Information](#pipeline-information)

### IRMA

<details markdown="1">
<summary>Output files</summary>

- `irma/<sample>`
  - `amended_consensus/`
    - Assembled gene segment consensus sequences: `*.fa`
  - `figures/`
    - Coverage and variants plot for gene segment: `*-coverageDiagram.pdf`
    - Heuristics graph for gene segment: `*-heuristics.pdf`
    - Gene segment variant phasing heatmap using *experimental enrichment* distances: `*-EXPENRD.pdf`
    - Gene segment variant phasing heatmap using *modified Jaccard* distances: `*-JACCARD.pdf`
    - Gene segment variant phasing heatmap using *mutual association* distances: `*-MUTUALD.pdf`
    - Gene segment variant phasing heatmap using *normalized joint probability* distances: `*-NJOINTP.pdf`
    - Read filtering, QC and assembly info plots: `READ_PERCENTAGES.pdf`
  - `intermediate/`
    - Intermediate analysis output files for each step in IRMA assembly.
  - `logs/`
    - Counts and scores for assembly, QC and read mapping: `*_log.txt`
    - Configuration file for IRMA analysis: `FLU-*.sh`
    - Table of IRMA execution parameters: `run_info.txt`
  - `matrices/`
    - Variant phasing matrices to construct heatmaps under `figures/`: `*.sqm`
  - `secondary/`
    - Secondary assemblies and unmatched reads.
  - `tables/`
    - Summary of gene segment paired-end merging stats, if applicable: `*-pairingStats.txt`
    - Summary coverage stats for assembly of gene segment: `*-coverage.txt`
    - Stats for every position and allele in assembly of gene segment: `*-allAlleles.txt`
    - Insertion variants called for gene segment: `*-insertions.txt`
    - Deletion variants called for gene segment: `*-deletions.txt`
    - SNP variants called for gene segment: `*-variants.txt`
    - Read counts at various stages of IRMA assembly process: `READ_COUNTS.txt`
  - Sorted BAM file for gene segment assembly: `*.bam`
  - Final assembled plurality consensus (no mixed basecalls) for gene segment: `*.fa`
  - IRMA variant call file for gene segment: `*.vcf`
- Concatenated "amended" IRMA consensus sequences for all gene segments assembled: `<sample>.irma.consensus.fasta`

</details>

[IRMA][] output is described in the [official IRMA output documentation](https://wonder.cdc.gov/amd/flu/irma/output.html).

The primary output from [IRMA][] are the consensus sequences for gene segments, which are used for H/N subtyping and performed blastn against influenza database to pull top match reference sequences for each segment of each sample.

### BLAST analysis

<details markdown="1">
<summary>Output files</summary>

- `blast/ncbi/blast_db/`
  - Nucleotide [BLAST] database of [NCBI Influenza DB][] and reference database (if provided option `--ref_db`): `influenza_db.*`
- `blast/ref_db/blast_db/`
  - Nucleotide [BLAST] database of the reference database (if provided option `--ref_db`) ref_fasta.fixed.*`
- `blast/blastn/irma`
  - Nucleotide [BLAST] tabular output files (`-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle"`) of sample IRMA assembled gene segments against the [NCBI Influenza DB][] and the reference database (if provided option `--ref_db`)
- `blast/blastn/consensus`
  - Nucleotide [BLAST] tabular output files (`-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle"`) of sample final consensus assembled gene segments against the [NCBI Influenza DB][] and the reference database (if provided option `--ref_db`)
- `blast/blastn/against_ref_db`
  - Nucleotide BLAST tabular output files (`-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle"`) of sample final consensus assembled gene segments against the reference database only (if provided option `--ref_db`)

</details>

Nucleotide [BLAST][] (`blastn`) is used to query [IRMA][] assembled gene segment sequences against [Influenza sequences from NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/) (and optionally, against user-specified sequences (`--ref_db`) to predict the H and N subtype of each sample if possible (i.e. if segments 4 (hemagglutinin) and/or 6 (neuraminidase) were assembled) and to determine the closest matching reference sequence for each segment for reference mapped assembly.

### Coverage Plots

<details markdown="1">
<summary>Output files</summary>

- `coverage_plots/<sample>/`
  - Coverage plot in linear and log scale: `*.pdf`

</details>

### Assembled Consensus Sequences

<details markdown="1">
<summary>Output files</summary>

- `consensus/bcftools/<sample>/`
  - Assembled consensus sequences for each segment: `*.bcftools.consensus.fasta`
- `consensus/bcftools/`
  - Concatenated consensus sequences for all segments assembled: `<sample>.consensus.fasta`
- `consensus/irma/`
  - Assembled consensus sequences for each segment: `<sample>.irma.consensus.fasta`

</details>

### Mismatch Report

<details markdown="1">
<summary>Output files</summary>

- `mistmacth_report/`
  - `<sample>-blastn-report.xlsx`

</details>
The report contains 2 sheets:

- **1_Mismatch_Report**: Count number of mismatches in BLASTN report (see sheet 2) against each reference sequences in reference database
- **2_Blastn_Results**: Nucleotide BLASTN report of sample final consensus against reference database

### Reference Sequences

<details markdown="1">
<summary>Output files</summary>

- `<sample>/`
  - Top reference sequences for all segments: `*.reference.fasta`
  - List of top reference ID pulled from influenza database: `*.topsegments.csv`

</details>

### Segments Mapping

<details markdown="1">
<summary>Output files</summary>

- `mapping/<sample>/`
  - The results of segments mapping using minimap2: `*.bam`, `*.bai`, `*.depths.tsv`, `*.flagstat`, `*.idxstats`, `*.stats`

</details>

### Variant Calling

<details markdown="1">
<summary>Output files</summary>

- `variants/<sample>/`
  - Filter Frameshift VCF: `*.filt_frameshift.vcf`
  - BCF Tools stats: `*.bcf_tools.stats.txt`
  - Clair3, Medaka, or Freebayes output directory

</details>

### H/N Subtyping

<details markdown="1">
<summary>Output files</summary>

- H/N subtyping Excel report: `nf-flu-subtyping-report.xlsx`

</details>

A H/N subtyping Excel report is generated from all [BLAST analysis](#blast-analysis) results for all samples and final assembled gene segments. The H and N subtypes are based on the proportion of high-quality BLAST matches that support the subtype prediction, that is, the top BLAST match for the HA and NA segments does not determine the subtype since the metadata for the top match could be incorrectly entered into NCBI.

More specifically, H and N subtype is predicted based on determining what the subtype is using the BLAST analysis results starting at a % identity threshold of 99% and decrementing by 1% until a subtype or the minimum % identity is reached (default: 85%). At least 3 hits are required to determine a subtype at a particular threshold. If no subtype is determined, the subtype is set to "N/A".

The subtyping report spreadsheet contains four sheets:

- **1_Subtype Predictions**: H/N subtype prediction results for each sample along with top matching Influenza DB segment for the H and N segments
- **2_Top Segment Matches**: Top 3 Influenza DB sequence matches for each segment of each sample along with BLASTN hit values and reference sequence metadata.
- **3_H Segment Results**: Top H subtype prediction, BLASTN results and top matching sequence metadata for each sample.
- **4_N Segment Results**: Top N subtype prediction, BLASTN results and top matching sequence metadata for each sample.

#### Sheet: 1_Subtype Predictions

This sheet contains the H/N subtype prediction results for each sample along with top matching Influenza DB segment for the H and N segments

| Field | Description                                                                                                                                                                                                                                                                                               | Example |
|-------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| Sample | Sample name                                                                                                                                                                                                                                                                                               | ERR3338653 |
| Subtype Prediction | H/N subtype prediction based on [BLAST][] against the Influenza DB. If a type could not be assigned to either H or N segment or both, then the subtype prediction will be missing the H or N value or if both the H and N cannot be assigned then the subtype prediction will be null or an empty cell value | H1N1 |
| H: top match accession | NCBI accession of top matching Influenza sequence for the H segment                                                                                                                                                                                                                                       | CY147779 |
| H: type prediction | H subtype prediction number. Value is a number.                                                                                                                                                                                                                                                           | 1 |
| H: top match virus name | Top matching sequence virus name                                                                                                                                                                                                                                                                          | Influenza A virus (A/Mexico/24036/2009(H1N1)) |
| H: NCBI Influenza DB subtype match proportion | Proportion of BLAST matches that support the H subtype prediction. This value is a decimal number where 1.0 indicates 100% of matches support the subtype prediction.                                                                                                                                     | 0.9980057896 |
| N: top match accession | NCBI accession of top matching Influenza DB sequence for the N segment                                                                                                                                                                                                                                    | MN371610 |
| N: type prediction | N subtype prediction. Value is a number.                                                                                                                                                                                                                                                                  | 1 |
| N: top match virus name  | Top matching sequence virus name                                                                                                                                                                                                                                                                          | Influenza A virus (A/California/04/2009) |
| N: NCBI Influenza DB subtype match proportion | Proportion of BLAST matches that support the N subtype prediction. This value is a decimal number where 1.0 indicates 100% of matches support the subtype prediction.                                                                                                                                     | 0.9993240503 |

#### Sheet: 2_Top Segment Matches

This sheet contains the top N Influenza DB sequence matches for each segment of each sample along with BLASTN hit values and reference sequence metadata.

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

### Annotation

Consensus sequences are annotated using [VADR][]. The output files are available in multiple formats including Feature Table, Genbank, nucleotide and amino acid FASTA and GFF.

<details markdown="1">
<summary>Output files</summary>

- `annotation/vadr/<sample>`/`annotation/vadr/irma/<sample>`/`annotation/vadr/bcftools/<sample>`
  - Each sample will have its own VADR annotation analysis output directory. Feature table output can be found in the `*.vadr.pass.tbl` files.
- `annotation/<sample>`/`annotation/vadr/irma/<sample>`/`annotation/vadr/bcftools/<sample>`
  - VADR Feature Table output is converted to Genbank, GFF and FASTA format for downstream analyses. FASTA files with nucleotide sequences of genetic features (CDS, mature peptide, signal peptide, etc) can be found in the `.ffn` files and amino acid sequences of genetic features can be found in the `.faa` files.
- `annotation/vadr-annotation-failed-sequences.txt`: list of sequences that failed VADR annotation
- `annotation/vadr-annotation-issues.txt`: table describing sequences that had issues with VADR annotation

</details>

### FluMut

Consensus sequences are analyzed using [FluMut] to identify mutations of interest, specifically for H5N1. The output files are available in multiple formats including Excel and TSV.

<details markdown="1">
<summary>Output files</summary>

- `flumut/`
  - `flumut.xlsm`: Excel file containing all the results from the FluMut analysis.
  - `flumut-markers.tsv`: Tab-separated file containing the list of detected markers by FluMut analysis.
  - `flumut-mutations.tsv`: Tab-separated file containing the list of amino acids present in the positions of mutations of interest for each sample by FluMut analysis.
  - `seqs-for-flumut.fasta`: FASTA file containing the sequences of the samples that were analyzed by FluMut.

</details>

### Genin2

[Genin2] genotypes consensus sequences based on information from clade 2.3.4.4b H5Nx viruses collected in Europe since October 2020.

<details markdown="1">
<summary>Output files</summary>

- `genin2/`
  - `genin2.tsv`: Tab-separated file containing the Genin2 genotyping results.

</details>

### GenoFLU

[GenoFlu](https://github.com/USDA-VS/GenoFLU/) "was developed to classify HPAI H5N1 goose/Guangdong clade 2.3.4.4b viruses detected in North American flyways. This tool considers all eight gene segments and can classify clade 2.3.4.4b viruses that have reassorted with North American low pathogenic viruses. The GenoFLU tool was developed for North America utilizing references detected primarily within the United States. The A1 GenoFLU genotype corresponds to the European National Reference Laboratory (EURL) genotype “C”, which is Eurasian wigeon/Netherlands-like virus that was predominant at the time the A1 virus was initially identified in Newfoundland."

<details markdown="1">
<summary>Output files</summary>

- `genoflu/`
  - `<sample>.genoflu.tsv`: Tab-separated file containing the results of the GenoFLU analysis for each sample.
  - `<sample>.genoflu.xlsx`: Excel file containing the results of the GenoFLU analysis for each sample.

</details>

### Cleavage Site Prediction

HA cleavage site prediction is performed using `bin/cleavage_site.py` with the [VADR][] annotation sequences.
The script will also classify the pathogenicity of the cleavage site based on the predicted cleavage site sequence.

<details markdown="1">
<summary>Output files</summary>

- `cleavage/`
  - `<sample>.cleavage.tsv`: Tab-separated file containing the results of the cleavage site prediction analysis for each sample.

</details>

### Nextclade

[Nextclade][] performs clade assignment, mutation calling and sequence quality checks.
In nf-flu, Nextclade analysis of assembled Influenza sequences is performed against 30 Nextclade datasets for different subtypes and lineages of Influenza A and B virus.

The specific Nextclade datasets and optionally versions (tags) can be configured with a headerless CSV file. Nextclade results are aggregrated across samples and datasets and filtered for positive results into a single Nextclade TSV (tab-separated values) report with additional fields capturing sample, dataset name and dataset version/tag information as well as Nextclade and Nextclade dataset specific results.

<details markdown="1">
<summary>Output files</summary>

- `nextclade/`
  - `nextclade.tsv`: Tab-separated file containing the results of Nextclade analysis across all assembled sample sequences and Nextclade datasets filtered for positive results.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.fixed.csv`.
  - Documentation for interpretation of results in HTML format: `results_description.html`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[FluMut]: https://github.com/izsvenezie-virology/FluMut
[Genin2]: https://github.com/izsvenezie-virology/genin2
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[NCBI Influenza DB]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[Nextclade]: https://clades.nextstrain.org/
[VADR]: https://github.com/ncbi/vadr
