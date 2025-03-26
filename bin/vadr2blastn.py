#!/usr/bin/env python

"""
Parsing of BLASTN output into a subtyping report expects that the FASTA headers
contain the sample name and the segment index. This script will parse the GenBank
and MDL (.mdl) files output by VADR and generate a FASTA file for BLASTN.


Example Influenza B virus VADR .mdl file contents:

#                                     num   num   num
#idx  model     group      subgroup  seqs  pass  fail
#---  --------  ---------  --------  ----  ----  ----
1     AF387493  fluB-seg4  -            1     1     0
2     AY191501  fluB-seg6  -            1     1     0
3     AY504599  fluB-seg2  -            1     1     0
4     AY504605  fluB-seg7  -            1     1     0
5     AY504614  fluB-seg8  -            1     0     1
6     EF626631  fluB-seg5  -            1     1     0
7     EF626633  fluB-seg3  -            1     1     0
8     EF626642  fluB-seg1  -            1     1     0
#---  --------  ---------  --------  ----  ----  ----
-     *all*     -          -            8     7     1
-     *none*    -          -            0     0     0
#---  --------  ---------  --------  ----  ----  ----


Example Influenza A virus VADR .mdl file contents:

#                                     num   num   num
#idx  model     group      subgroup  seqs  pass  fail
#---  --------  ---------  --------  ----  ----  ----
1     CY000449  fluA-seg4  H1           1     1     0
2     CY002009  fluA-seg7  -            1     1     0
3     CY002079  fluA-seg1  -            1     1     0
4     CY002284  fluA-seg8  -            1     1     0
5     CY002538  fluA-seg6  N1           1     1     0
6     CY003645  fluA-seg3  -            1     1     0
7     CY003646  fluA-seg2  -            1     1     0
8     CY006079  fluA-seg5  -            1     1     0
#---  --------  ---------  --------  ----  ----  ----
-     *all*     -          -            8     8     0
-     *none*    -          -            0     0     0
#---  --------  ---------  --------  ----  ----  ----

"""

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

FLU_GENE_SEGMENT_IDX = {
    "A": {
        "PB2": 1,
        "PB1": 2,
        "PA": 3,
        "HA": 4,
        "NP": 5,
        "NA": 6,
        "MP": 7,
        "NS": 8,
    },
    "B": {
        "PB1": 1,
        "PB2": 2,
        "PA": 3,
        "HA": 4,
        "NP": 5,
        "NA": 6,
        "MP": 7,
        "NS": 8,
    },
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


def get_vadr_mdl_flu_a_or_b(mdl_path: Path) -> str:
    """Determine if the sample is Influenza A or B virus based on VADR MDL file

    This file should contain lines containing "fluA-seg\\d" or "fluB-seg\\d" 
    to determine the species of influenza virus.
    """
    with open(mdl_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "fluA-seg" in line:
                return "A"
            if "fluB-seg" in line:
                return "B"
    return ''


def parse_influenza_gbk(
        sample_name: str,
        gbk_path: Path,
        flu_a_or_b: str
) -> Iterator[SeqRecord]:
    if flu_a_or_b not in ("A", "B"):
        raise ValueError(f"Invalid value for 'flu_a_or_b': {flu_a_or_b}")

    if sample_name == '':
        sample_name = gbk_path.stem

    for seq_record in SeqIO.parse(gbk_path, "genbank"):
        segment = get_segment(seq_record)
        if not segment:
            logger.warning(f"Could not determine segment for {seq_record.id}")
            continue
        segment_idx = FLU_GENE_SEGMENT_IDX[flu_a_or_b].get(segment)
        if not segment_idx:
            logger.warning(f"Could not determine segment index for {segment} of '{gbk_path}'")
            continue
        seq_record.id = f"{sample_name}_{segment_idx}"
        yield seq_record


def version_callback(value: bool) -> None:
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()


@app.command()
def main(
        sample_name: str,
        gbk_path: Path,
        vadr_outdir: Path,
        version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True, is_flag=True),
):
    """VADR GenBank and mdl to BLASTN FASTA

    VADR GenBank and MDL table output will be used to generate a FASTA file for 
    BLASTN where the segment index is appended to the sample name. The sample
    name is derived from the GenBank file name.
    """
    init_logging(verbose=True)
    mdl_path = next(vadr_outdir.glob("*.vadr.mdl"), None)
    if not mdl_path:
        logger.error(f"Could not find VADR .mdl file in '{vadr_outdir}'")
        raise typer.Exit(code=1)
    flu_a_or_b = get_vadr_mdl_flu_a_or_b(mdl_path)
    if not flu_a_or_b:
        logger.warning(f"Could not determine influenza A or B for {mdl_path}")
        return
    for seq_record in parse_influenza_gbk(sample_name, gbk_path, flu_a_or_b):
        print(f'>{seq_record.id}')
        print(str(seq_record.seq))


if __name__ == "__main__":
    app()
