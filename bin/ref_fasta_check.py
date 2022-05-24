#!/usr/bin/env python

from pathlib import Path
from Bio import SeqIO
import typer
import re
from rich.logging import RichHandler
import logging

def main(input_fasta: Path, output_fasta: Path):
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    logging.info(
        f'input_fasta="{input_fasta}" output_correct_fasta="{output_fasta}"'
    )
    ref_fasta = SeqIO.parse(open(input_fasta), 'fasta')
    with open(output_fasta, 'w') as outfile:
        for record in ref_fasta:
            seqid, sequence = record.id.strip(), record.seq
            seq_record_id = re.sub(r"[()\"#/@;:<>{}`+=~|!?,]", "_", seqid)
            outfile.write(f'>{seq_record_id}\n{sequence}\n')


if __name__ == "__main__":
    typer.run(main)
