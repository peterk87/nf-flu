# About `FB701709.1-with-end-snps.badread.bam`

Simulated Nanopore reads with [Badread](https://github.com/rrwick/Badread) from modified reference sequence for FB701709.1 where 2 SNP were introduced within the first and last 16bp (A5C, C2337A):

```
$ badread simulate --reference ../fasta/FB701709.1-with-end-snps.fasta --quantity 30x | pigz > FB701709.1-with-end-snps.badread.fastq.gz

Badread v0.4.1
long read simulation

Loading reference from ../fasta/FB701709.1-with-end-snps.fasta
  1 contig:
    FB701709.1: 2,341 bp, linear, 1.00x depth

Generating fragment lengths from a gamma distribution:
  mean  =  15000 bp      parameters:
  stdev =  13000 bp        k (shape)     = 1.3314e+00
  N50   =  22622 bp        theta (scale) = 1.1267e+04
  │  ▖▌▌▌▌▌▖▖▖
  │ ▖▌▌▌▌▌▌▌▌▌▌▌▖▖
f │ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖
r │▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖
a │▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖
g │▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖
s │▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖▖▖
  │▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖▖▖▖▖▖▖▖▖▖▖▖▖▖▖▖
  ├─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
  │               ▖▖▖▌▌▌▌▌▌▌▌▌▖▖▖▖▖
  │           ▖▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖
b │         ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖
a │       ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖
s │     ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖
e │    ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖▖
s │  ▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖▖▖▖▖
  │ ▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▖▖▖▖▖▖▖
  ├─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
  0         6800      13600     20400     27200     34000     40800     47600     54400     61200     68000

Generating read identities from a beta distribution:
  mean  =  95%      shape parameters:
  max   =  99%        alpha = 5.7384e+01
  stdev = 2.5%        beta  = 2.4162e+00
  │                                                                                           ▖▌▌▖
  │                                                                                          ▌▌▌▌▌
  │                                                                                        ▖▌▌▌▌▌▌▌
  │                                                                                       ▖▌▌▌▌▌▌▌▌
  │                                                                                      ▌▌▌▌▌▌▌▌▌▌▌
  │                                                                                    ▖▌▌▌▌▌▌▌▌▌▌▌▌
  │                                                                                 ▖▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌
  │                                                                           ▖▖▖▖▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌▌
  ├─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┬─────────┐
  50        55        60        65        70        75        80        85        90        95        100

Loading error model from /home/user/miniforge3/envs/badread/lib/python3.12/site-packages/badread/error_models/nanopore2023.gz
  done: loaded error distributions for 16384 7-mers

Loading qscore model from /home/user/miniforge3/envs/badread/lib/python3.12/site-packages/badread/qscore_models/nanopore2023.gz
  done: loaded qscore distributions for 10000 alignments

Read glitches:
  rate (mean distance between glitches) = 10000
  size (mean length of random sequence) =    25
  skip (mean sequence lost per glitch)  =    25

Start adapter:
  seq: AATGTACTTCGTTCAGTTACGTATTGCT
  rate:   90.0%
  amount: 60.0%

End adapter:
  seq: GCAATACGTAACTGAACGAAGT
  rate:   50.0%
  amount: 20.0%

Other problems:
  chimera join rate: 1%
  junk read rate:    1%
  random read rate:  1%

Target read set size: 70,230 bp

Simulating: 35 reads  72,294 bp  100.0%
```

Reads were aligned to original FB701709.1 sequence with Minimap2:

```bash
$ minimap2 -t 4 -ax map-ont ../fasta/FB701709.1.fasta FB701709.1-with-end-snps.badread.fastq.gz | samtools sort | samtools view -F4 -b > FB701709.1-with-end-snps.badread.bam
```

Testing that variants within the 16bp of each end are called by Clair3 with the new `--enable_variant_calling_at_sequence_head_and_tail` option:

```
$ run_clair3.sh \
  --haploid_sensitive \
  --enable_long_indel \
  --keep_iupac_bases \
  --include_all_ctgs \
  --var_pct_phasing=1 \
  --var_pct_full=1 \
  --ref_pct_full=1 \
  --platform="ont" \
  --enable_variant_calling_at_sequence_head_and_tail \
  --bam_fn=FB701709.1-with-end-snps.badread.bam \
  --ref_fn=../fasta/FB701709.1.fasta \
  --model_path="$(dirname $(which run_clair3.sh))/models/r941_prom_sup_g5014" \
  --output clair3-out
```

The `merge_output.vcf.gz` should contain the 2 artificially introduced SNPs:

```
##fileformat=VCFv4.2
##source=Clair3
##clair3_version=1.0.11
##cmdline=/home/user/miniforge3/envs/clair3/bin/run_clair3.sh --haploid_sensitive --enable_long_indel --keep_iupac_bases --include_all_ctgs --var_pct_phasing=1 --var_pct_full=1 --ref_pct_full=1 --platform=ont --enable_variant_calling_at_sequence_head_and_tail --bam_fn=FB701709.1-with-end-snps.badread.bam --ref_fn=../fasta/FB701709.1.fasta --model_path=/home/user/miniforge3/envs/clair3/bin/models/r941_prom_sup_g5014 --threads=4 --output clair3-out
##reference=/home/user/repos/nf-flu/tests/data/reads/../fasta/FB701709.1.fasta
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality variant">
##FILTER=<ID=RefCall,Description="Reference call">
##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">
##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">
##contig=<ID=FB701709.1,length=2341>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
FB701709.1      5       .       A       C       22.77   PASS    P       GT:GQ:DP:AD:AF  1:22:25:0,25:1.0000
FB701709.1      2337    .       C       A       15.91   PASS    F       GT:GQ:DP:AD:AF  1:15:23:0,23:1.0000
```
