#!/usr/bin/env python

import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd
import typer
from Bio.SeqIO.FastaIO import SimpleFastaParser
from rich.logging import RichHandler

VERSION = "2025.01.1"

logger = logging.getLogger(__name__)
app = typer.Typer(rich_markup_mode="markdown")


def init_logging(verbose: bool) -> None:
    from rich.traceback import install
    install(show_locals=True, width=120, word_wrap=True)

    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG if verbose else logging.INFO,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def get_ha_seqs_and_features(input_fasta: Path) -> dict[str, dict[str, tuple[str, str]]]:
    out = defaultdict(dict)
    with open(input_fasta) as fh:
        for header, seq in SimpleFastaParser(fh):
            annotation = ""
            if header.endswith("HA1"):
                annotation = "HA1"
            elif header.endswith("HA2"):
                annotation = "HA2"
            elif header.endswith("hemagglutinin"):
                annotation = "hemagglutinin"
            elif header.endswith("cleavage site"):
                annotation = "cleavage site"
            if not annotation:
                continue
            header_first_part, *_ = header.split(" ")
            seqid = header_first_part.split("|")[0]
            seqid = seqid.replace('_segment4_HA', '')
            seq = seq.upper()
            out[seqid][annotation] = (header, seq)
    return out


def find_cleavage_sites(sample_ha1_ha2_seqs: dict[str, dict[str, tuple[str, str]]]) -> list[dict] | None:
    # If no pairs were passed, log and exit
    if not sample_ha1_ha2_seqs:
        logger.error("No cleavage sequences found. Please check the input file or extraction logic.")
        return None

    out = []
    for sample, ha_dict in sample_ha1_ha2_seqs.items():
        cleavage_site_header, cleavage_site_seq = ha_dict.get("cleavage site", (None, None))
        
        basic_residue_count = 0
        if cleavage_site_seq:
            basic_residue_count, pathogenicity, residue_type = count_ha1_basic_residues(cleavage_site_seq)
            logger.info(
                f"{sample}: cleavage sequence '{cleavage_site_seq}': {residue_type=}; {basic_residue_count=}; {pathogenicity=}")
        else:
            logger.info(f"{sample}: No cleavage site found")
            residue_type = "N/A"
            pathogenicity = "N/A"

        cleavage_site_search_result = {
            'sample': sample,
            'cleavage_seq': cleavage_site_seq or "No cleavage site found",
            'residue_type': residue_type,
            'basic_residue_count': basic_residue_count,
            'pathogenicity': pathogenicity,
            'cleavage_site_header': cleavage_site_header or "N/A",
        }
        out.append(cleavage_site_search_result)
        logger.debug(f"{cleavage_site_search_result=}")
    return out


def count_ha1_basic_residues(ha_cleavage_site_seq: str) -> tuple[int, str, str]:
    """Count basic residues in HA1 motif and classify based on count.

    >>> count_ha1_basic_residues("PQLLLLLLLLLRGLF")
    (1, 'LPAI', 'Monobasic')
    >>> count_ha1_basic_residues("PEIPKRRRRGIF")
    (5, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PEIPKRRGLF")
    (3, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQIESRGLF")
    (1, 'LPAI', 'Monobasic')
    >>> count_ha1_basic_residues("PQGERRRKKRGLF")
    (6, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQERREKRGLF")
    (4, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQERERGLF")
    (2, 'LP/HP', 'Multibasic')
    >>> count_ha1_basic_residues("PQERREERGLF")
    (1, 'LPAI', 'Monobasic')
    >>> count_ha1_basic_residues("R")
    (0, 'N/A', 'N/A')
    >>> count_ha1_basic_residues("RK")
    (0, 'N/A', 'N/A')
    >>> count_ha1_basic_residues("RKK")
    (0, 'N/A', 'N/A')
    >>> count_ha1_basic_residues("RRKK")
    (0, 'N/A', 'N/A')
    >>> count_ha1_basic_residues("RRRKK")
    (1, 'LPAI', 'Monobasic')
    """
    if not ha_cleavage_site_seq:
        return 0, "N/A", "N/A"
    if len(ha_cleavage_site_seq) < 4:
        return 0, "N/A", "N/A"
    if ha_cleavage_site_seq[-4] not in "RK":
        return 0, "N/A", "N/A"
    # trim last 3 residues
    ha_cleavage_site_seq = ha_cleavage_site_seq[:-3]
    basic_residue_count = 0
    consecutive_non_basic = 0
    for i in range(len(ha_cleavage_site_seq) - 1, 0, -1):
        residue = ha_cleavage_site_seq[i]
        if residue in ("R", "K"):
            basic_residue_count += 1
            consecutive_non_basic = 0
        else:
            consecutive_non_basic += 1
            # if more than one non-basic residue is found before the last 3 residues, break
            # or if a non-basic residue is found after the last 3 residues, break
            if (
                    i > len(ha_cleavage_site_seq) - 4
                    and consecutive_non_basic > 1
                    or i <= len(ha_cleavage_site_seq) - 4
            ):
                break
    # Classify based on basic residue count
    if basic_residue_count == 1:
        residue_type = "Monobasic"
    elif basic_residue_count > 1:
        residue_type = "Multibasic"
    else:
        residue_type = "N/A"
    if basic_residue_count == 1:
        pathogenicity = "LPAI"
    elif 1 < basic_residue_count <= 2:
        pathogenicity = "LP/HP"
    elif basic_residue_count >= 3:
        pathogenicity = "HPAI"
    else:
        pathogenicity = "N/A"
    return basic_residue_count, pathogenicity, residue_type


def write_sites(output_file: Path, sites: list[dict]) -> None:
    """Process the sequences to detect cleavage sites and classify them."""
    cols = [
        ('sample', 'Sample'),
        ('cleavage_seq', 'Cleavage Sequence'),
        ('residue_type', 'Residue Type'),
        ('basic_residue_count', 'Basic Residue Count'),
        ('pathogenicity', 'Pathogenicity'),
        ('cleavage_site_header', 'Cleavage Site Sequence Header'),
    ]

    df = pd.DataFrame(sites)
    df.rename(columns=dict(cols), inplace=True)
    df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"Results written to: {output_file}")


def version_callback(value: bool):
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()


@app.command()
def main(
        input_fasta: Path = typer.Argument(
            ...,
            exists=True,
            help="Input Influenza amino acid sequences FASTA file. "
                 "Preferably the output from VADR annotation post-table2asn."
        ),
        output_tsv: Path = typer.Argument(Path("cleavage_site.tsv"), help="Output TSV file"),
        force: bool = typer.Option(False, "--force", "-f", help="Overwrite output file if it exists"),
        verbose: bool = typer.Option(False, "--verbose", "-v", is_flag=True, help="Enable verbose logging"),
        version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
):
    """# cleavage_site.py: classify Influenza A virus HA cleavage sites from VADR annotation amino acid sequences.

    The script reads a VADR annotation FASTA file amino acid sequences and classifies HA cleavage sites based on the
    number of basic residues in the motif, i.e. 'Monobasic', 'Multibasic', or 'N/A' and pathogenicity ('LPAI', 'LP/HP',
    'HPAI', or 'N/A').

    ## Output

    The results are written to a TSV file with the following columns:

    * **Sample**: Sample ID

    * **Cleavage Sequence**: Detected cleavage sequence (HA1 motif matching sequence + first 3 HA2 residues)

    * **Residue Type**: Classification of the basic residues in the HA1 motif

    * **Basic Residue Count**: Number of basic residues in the HA1 motif

    * **Pathogenicity**: Classification of the pathogenicity based on the basic residue count
    """
    init_logging(verbose)
    if output_tsv.exists() and not force:
        logger.error(f"Output file '{output_tsv}' exists. Use --force to overwrite.")
        raise typer.Abort()
    logger.info(f"{input_fasta=}")
    # Step 1: Extract HA1 and HA2 sequences
    ha_seqs = get_ha_seqs_and_features(input_fasta)
    logger.info(f"Read {len(ha_seqs)} HA sequences and features from '{input_fasta}'")
    logger.debug(f"{ha_seqs=}")
    # Step 2: Process cleavage site analysis and classification
    sites = find_cleavage_sites(ha_seqs)
    write_sites(output_tsv, sites)


if __name__ == "__main__":
    app()
