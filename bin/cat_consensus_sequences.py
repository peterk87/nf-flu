#!/usr/bin/env python

from Bio import SeqIO
import click


@click.command()
@click.option("-s", "--sample-name", default="", help="Sample Name.")
@click.option("-o1", "--output1-fasta", default="consensus.fasta", help="Consensus Fasta")
@click.option("-o2", "--output2-fasta", default="blastn.fasta", help="Consensus Fasta for Blastn")
@click.argument("consensus_sequences", nargs=-1)
def write_consensus(sample_name, output1_fasta, output2_fasta, consensus_sequences):
    sequence_list = []
    for infile in consensus_sequences:
        fasta_file = SeqIO.parse(open(infile), 'fasta')
        for record in fasta_file:
            seqid, sequence = record.id, record.seq
            *_, segment_number, segment_name = seqid.split('_')
            sequence_list.append([seqid, segment_number, segment_name, sequence])
    sequence_list.sort(key=lambda tup: tup[1])
    # Outfile for publishing to output dir with header format SampleName_segment1_PB2
    with open(output1_fasta, 'w') as outfile1:
        for item1 in sequence_list:
            outfile1.write(f'>{sample_name}_segment{item1[1]}_{item1[2]}\n{item1[3]}\n')
    # Outfile for blastn report script compatible
    with open(output2_fasta, 'w') as outfile2:
        for item2 in sequence_list:
            outfile2.write(f'>{sample_name}_{item2[1]}\n{item2[3]}\n')


if __name__ == "__main__":
    write_consensus()
