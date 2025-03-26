#!/usr/bin/env python

import logging
from pathlib import Path
from typing import Iterator

import typer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rich.logging import RichHandler

VERSION = "2025.03.1"

logger = logging.getLogger(__name__)
app = typer.Typer(rich_markup_mode="markdown")

GENE_TO_SEGMENT = {
    "PB2": "PB2",
    "PB1": "PB1",
    # segment 3
    "PA": "PA",
    # segment 4
    "HA": "HA",
    "HA1": "HA",
    "HA2": "HA",
    # segment 5
    "NP": "NP",
    # segment 6
    "NA": "NA",
    # segment 7
    "M1": "MP",
    "M2": "MP",
    # segment 8
    "NEP": "NS",
    "NS1": "NS",
    "NS2": "NS",
}


def init_logging(verbose: bool) -> None:
    from rich.traceback import install
    install(show_locals=True, width=120, word_wrap=True)

    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG if verbose else logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def get_segment(seq_record: SeqRecord) -> str:
    genes = [f.qualifiers["gene"][0] for f in seq_record.features if "gene" in f.qualifiers]
    if not genes:
        return ''
    segments = [GENE_TO_SEGMENT.get(gene, '') for gene in genes]
    # remove empty elements
    segments = [segment for segment in segments if segment]
    # check if all elements in segments are the same
    if len(set(segments)) == 1:
        return segments[0]
    # otherwise return the first non-empty element
    for segment in segments:
        if segment:
            return segment
    return ''


def parse_influenza_gbk(gbk_path: Path) -> Iterator[SeqRecord]:
    sample_name = gbk_path.stem
    for seq_record in SeqIO.parse(gbk_path, "genbank"):
        segment = get_segment(seq_record)
        seq_record.id = f"{sample_name}_{segment}"
        yield seq_record


def version_callback(value: bool) -> None:
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()


@app.command()
def main(
        gbks_dir_path: Path,
        version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True, is_flag=True),
):
    """VADR GenBank to FluMut FASTA

    Expecting VADR Genbank files in a directory where each sample has a single GenBank file. 
    The sample name is the GenBank file name without the extension and will be used within the
    output FASTA headers along with the VADR annotation determined segments.
    """
    init_logging(verbose=True)
    gbks_dir_path = Path(gbks_dir_path)
    if not gbks_dir_path.is_dir():
        logger.error(f"Invalid directory path: {gbks_dir_path}")
        raise typer.Exit(code=1)

    for gbk_path in gbks_dir_path.glob("*.gbk"):

        for seq_record in parse_influenza_gbk(gbk_path):
            print(f'>{seq_record.id}')
            print(str(seq_record.seq))


if __name__ == "__main__":
    app()
