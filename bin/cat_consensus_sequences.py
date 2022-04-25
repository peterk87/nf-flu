#!/usr/bin/env python

from pathlib import Path
from Bio import SeqIO
import typer
import click
import re
from typing import Dict, List, Optional, Tuple
from rich.logging import RichHandler
import logging

@click.command()
@click.option("-s", "--sample-name", default="", help="Sample Name.")
@click.option("-o", "--output-fasta", default="consensus.fasta", help="Consensus Fasta")
@click.argument("consensus_sequences", nargs=-1)
def write_consensus(sample_name, output_fasta, consensus_sequences):

    sequence_list = []
    for infile in consensus_sequences:
        fasta_file = SeqIO.parse(open(infile), 'fasta')
        for record in fasta_file:
            seqid, sequence = record.id.strip(), record.seq
            segment_number = seqid.split("_")[-2]
            sequence_list.append([seqid,segment_number,sequence])
    sequence_list.sort(key=lambda tup: tup[1])
    with open(output_fasta, 'w') as outfile:
        for record in sequence_list:
            outfile.write(f'>{sample_name+"_"+record[1]}\n{record[2]}\n')


if __name__ == "__main__":
    write_consensus()
