# Update Influenza sequences for nf-flu

Download the latest `nt_viruses` BLAST DB from NCBI using `update_blastdb.pl` and fetch the metadata table for all viruses in NCBI:

```bash
mkdir -p blastn-nt_viruses && cd blastn-nt_viruses
update_blastdb.pl --verbose --decompress nt_viruses

export DATE=$(date -I)
export projdir = "~/projects/$DATE-influenza-db-update" 
mkdir $projdir && cd $projdir
curl -SLk https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv.gz | gzip -cd > "$DATE-all-viruses-nucl-metadata.csv"
```

Using `dnaio`, `polars` and `blastdbcmd`, get Influenza sequences and their metadata for output  

```python
import polars as pl
from pathlib import Path

dfl = pl.scan_csv(Path("2025-04-04-all-viruses-nucl-metadata.csv"))

cols = [
    '#Accession',
    'Release_Date',
    'Genus',
    'Length',
    'Genotype',
    'Segment',
    'Geo_Location',
    'Host', 
    'Isolation_Source',
    'Collection_Date',
    'GenBank_Title'
]

flu_family = ['Orthomyxoviridae']
flu_genera = ['Alphainfluenzavirus', 'Betainfluenzavirus', 'Deltainfluenzavirus', 'Gammainfluenzavirus']

df = dfl.filter(
    pl.col('Family').is_in(flu_family) & 
    pl.col('Genus').is_in(flu_genera) &
    (pl.col('Length') >= 500)  
).collect()

df_out = df[cols]

df_out.write_csv('2025-04-04-influenza.csv')

with open('2025-04-04-influenza-accessions.txt', 'w') as fout:
    for acc in df['#Accession']:
        fout.write(f"{acc}\n")
```

The BLAST `blastdbcmd` may skip several thousand sequences, however, most of those sequences will be patent sequences and therefore not present within the `nt_viruses` BLAST DB.

```bash
blastdbcmd -entry_batch 2025-04-04-influenza-accessions.txt -db ~/blast-nt_viruses/nt_viruses -out 2025-04-04-influenza.fasta
```

This command will output the following message for several thousand likely patent sequences:

```text
Error: [blastdbcmd] Skipped PV440337.1
```

Select and output metadata for sequences that were successfully retrieved from the `nt_viruses` BLAST DB and ensure that only unique and non-redundant sequences are output for use with nf-flu:

```python
import dnaio
from pathlib import Path
import polars as pl
from tqdm import tqdm

df_out = pl.read_csv('2025-04-04-influenza.csv')

from collections import defaultdict
seq_hashes = defaultdict(list)
fasta_ids = []
with dnaio.open(Path('2025-04-04-influenza.fasta')) as reader:
    for record in reader:
        fasta_ids.append(record.id)
        seqhash = hash(str(record.sequence).upper())
        seq_hashes[seqhash].append(record.id)

seqids = set()
for xs in seq_hashes.values():
    seqid = xs[0]
    if seqid in seqids:
        print(len(xs))
    seqids.add(seqid)


with open('nr-influenza-seq-acc.txt', 'w') as fout:
    for seqid in seqids:
        fout.write(f'{seqid}\n')

df_subset = df_out.filter(pl.col('#Accession').is_in(seqids))

subset_acc_set = set(df_subset['#Accession'])

seqhashes = set()
seqids_set = set()
with dnaio.open(Path('2025-04-04-influenza.clean.fasta')) as reader, dnaio.open('2025-04-04-influenza.nr.fasta', mode='w', fileformat='fasta') as writer:
    for record in tqdm(reader):
        seqid = record.id
        seq = str(record.sequence).upper()
        seqhash = hash(seq)
        if seqhash in seqhashes:
            continue
        if seqid in seqids_set:
            continue
        if seqid not in subset_acc_set:
            continue
        seqhashes.add(seqhash)
        seqids_set.add(seqid)
        new_rec = dnaio.Sequence(name=seqid, sequence=seq)
        writer.write(new_rec)

df_subset.write_csv('2024-04-04-influenza.csv')
```

Zstandard compression of the FASTA and CSV files:

```bash
zstandard 2024-04-04-influenza.nr.fasta
zstandard 2024-04-04-influenza.csv
```
