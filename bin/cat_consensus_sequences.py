#!/usr/bin/env python

from Bio import SeqIO
import click


@click.command()
@click.option("-s", "--sample-name", default="", help="Sample Name.")
@click.option("-o1", "--output1-fasta", default="consensus.fasta", help="Consensus Fasta")
@click.option("-o2", "--output2-fasta", default="blastn.fasta", help="Consensus Fasta for Blastn")
@click.argument("fastas", nargs=-1)
def write_consensus(sample_name, output1_fasta, output2_fasta, fastas):
    seqs = []
    for fasta_path in fastas:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            *_, segment_number, segment_name = record.id.split('_')
            seqs.append([segment_number, segment_name, str(record.seq)])
    seqs.sort(key=lambda tup: tup[0])
    # Outfile for publishing to output dir with header format SampleName_segment1_PB2
    with open(output1_fasta, 'w') as fout:
        for segment_number, segment_name, seq in seqs:
            fout.write(f'>{sample_name}_segment{segment_number}_{segment_name}\n{seq}\n')
    # Outfile for blastn report script compatible
    with open(output2_fasta, 'w') as fout:
        for segment_number, _, seq in seqs:
            fout.write(f'>{sample_name}_{segment_number}\n{seq}\n')


if __name__ == "__main__":
    write_consensus()
