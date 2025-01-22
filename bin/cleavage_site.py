#!/usr/bin/env python

import logging
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import typer
from Bio.SeqIO.FastaIO import SimpleFastaParser
from rich.logging import RichHandler

VERSION = "2025.01.1"
HA1_REGEX = re.compile(r"""
    P            # Proline
    (Q|L|E)      # Glutamine, Leucine, or Glutamic Acid
    [A-Z]{1,10}  # Any character repeated 1-10 times
    (R|K)        # Arginine or Lysine
    (.*)         # Any character
""", re.VERBOSE)

HA1_REPLACE_REGEX = re.compile(r"""
    P     # Proline
    .     # Any character
    (R|K) # Arginine or Lysine
""", re.VERBOSE)

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


def get_ha1_ha2_seqs(input_fasta: Path) -> dict[str, dict[str, tuple[str, str]]]:
    out = defaultdict(dict)
    with open(input_fasta) as fh:
        for header, seq in SimpleFastaParser(fh):
            if not (header.endswith("HA1") or header.endswith("HA2")):
                continue
            header_first_part, *_, header_last_part = header.split(" ")
            seqid = header_first_part.split("|")[0]
            seqid = seqid.replace('_segment4_HA', '')
            seq = seq.upper()
            out[seqid][header_last_part] = (header, seq)
    return out


def find_cleavage_sites(sample_ha1_ha2_seqs: dict[str, dict[str, tuple[str, str]]]) -> list[dict] | None:
    # If no pairs were passed, log and exit
    if not sample_ha1_ha2_seqs:
        logger.error("No cleavage sequences found. Please check the input file or extraction logic.")
        return None

    out = []
    for sample, ha1_ha2_dict in sample_ha1_ha2_seqs.items():
        ha1_header, ha1_sequence = ha1_ha2_dict.get("HA1", (None, None))
        ha2_header, ha2_sequence = ha1_ha2_dict.get("HA2", (None, None))
        if not (ha1_sequence and ha2_sequence):
            logger.warning(f"HA1 or HA2 sequence missing for sample {sample}. Skipping.")
            continue
        logger.debug(f"{sample=} {ha1_header=} {ha2_header=}")
        logger.debug(f"{ha1_sequence=}")
        logger.debug(f"{ha2_sequence=}")
        # look at last 15 residues of HA1
        ha1_subseq = ha1_sequence[-15:]
        # first 3 residues of HA2
        ha2_subseq = ha2_sequence[:3]
        ha1_match = HA1_REGEX.search(ha1_subseq)
        basic_residue_count = 0
        if ha1_match:
            ha1_motif_match = ha1_match[0]
            cleavage_sequence = f"{ha1_motif_match}|{ha2_subseq}"
            # Remove 'P.(R|K)' pattern from HA1 for residue counting
            ha1_motif_match = HA1_REPLACE_REGEX.sub("", ha1_motif_match)
            basic_residue_count, pathogenicity, residue_type = count_ha1_basic_residues(ha1_motif_match)
            logger.info(
                f"{sample}: cleavage sequence found '{cleavage_sequence}' ({residue_type=}; {basic_residue_count=}; {pathogenicity=})")
        else:
            logger.info(f"{sample}: No cleavage site detected in HA1 sequence.")
            cleavage_sequence = "No cleavage site found"
            residue_type = "No match"
            pathogenicity = "check seq"

        cleavage_site_search_result = {
            'sample': sample,
            'cleavage_seq': cleavage_sequence,
            'residue_type': residue_type,
            'basic_residue_count': basic_residue_count,
            'pathogenicity': pathogenicity,
            'ha1_subseq': ha1_subseq,
            'ha2_subseq': ha2_subseq,
            'ha1_header': ha1_header,
            'ha1_sequence': ha1_sequence,
            'ha2_header': ha2_header,
            'ha2_sequence': ha2_sequence,
        }
        out.append(cleavage_site_search_result)
        logger.debug(f"{cleavage_site_search_result=}")
    return out


def count_ha1_basic_residues(ha1_motif_match: str) -> tuple[int, str, str]:
    """Count basic residues in HA1 motif and classify based on count.

    >>> count_ha1_basic_residues("PQLLLLLLLLLR")
    (1, 'LPAI', 'Monobasic')
    >>> count_ha1_basic_residues("PEIPKRRRR")
    (5, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PEIPKRR")
    (3, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQIESR")
    (1, 'LPAI', 'Monobasic')
    >>> count_ha1_basic_residues("PQGERRRKKR")
    (6, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQERREKR")
    (4, 'HPAI', 'Multibasic')
    >>> count_ha1_basic_residues("PQERER")
    (2, 'LP/HP', 'Multibasic')
    >>> count_ha1_basic_residues("PQERREER")
    (1, 'LPAI', 'Monobasic')
    """
    basic_residue_count = 0
    residue_type = ""
    # Count basic residues in HA1 sequence
    if ha1_motif_match[-1] not in "RK":
        logger.info(
            f"Last position ({len(ha1_motif_match)}:{ha1_motif_match[-1]}) is not basic! Labeling as 'check seq'.")
        residue_type = "check seq"
        basic_residue_count = 0
    else:
        consecutive_non_basic = 0
        for i in range(len(ha1_motif_match) - 1, 0, -1):
            residue = ha1_motif_match[i]
            if residue in ("R", "K"):
                basic_residue_count += 1
                consecutive_non_basic = 0
            else:
                consecutive_non_basic += 1
                # if more than one non-basic residue is found before the last 3 residues, break
                # or if a non-basic residue is found after the last 3 residues, break
                if (
                        i > len(ha1_motif_match) - 4
                        and consecutive_non_basic > 1
                        or i <= len(ha1_motif_match) - 4
                ):
                    break
        # Classify based on basic residue count
        if basic_residue_count == 1:
            residue_type = "Monobasic"
        elif basic_residue_count > 1:
            residue_type = "Multibasic"
        else:
            residue_type = "No match"
    if basic_residue_count == 1:
        pathogenicity = "LPAI"
    elif basic_residue_count == 2:
        pathogenicity = "LP/HP"
    elif basic_residue_count >= 3:
        pathogenicity = "HPAI"
    else:
        pathogenicity = "check seq"
    return basic_residue_count, pathogenicity, residue_type


def write_sites(output_file: Path, sites: list[dict]) -> None:
    """Process the sequences to detect cleavage sites and classify them."""
    cols = [
        ('sample', 'Sample'),
        ('cleavage_seq', 'Cleavage Sequence'),
        ('residue_type', 'Residue Type'),
        ('basic_residue_count', 'Basic Residue Count'),
        ('pathogenicity', 'Pathogenicity'),
        ('ha1_subseq', 'HA1 Subsequence (Last 15 residues)'),
        ('ha2_subseq', 'HA2 Subsequence (First 3 residues)'),
        ('ha1_header', 'HA1 Sequence Header'),
        ('ha1_sequence', 'HA1 Amino Acid Sequence'),
        ('ha2_header', 'HA2 Sequence Header'),
        ('ha2_sequence', 'HA2 Amino Acid Sequence'),
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
    """# cleavage_site.py: Find and classify Influenza A virus HA cleavage sites from amino acid sequences.

    The script reads a FASTA file with HA1 and HA2 amino acid sequences and searches for the HA1 cleavage site motif.
    The last 15 residues of HA1 are searched for the motif (regex: `P(Q|L|E)[A-Z]{1,10}(R|K)(.*)`). The number of basic
    residues in the motif is counted and classified as 'Monobasic', 'Multibasic', or 'No match'. The pathogenicity is
    classified as 'LPAI', 'LP/HP', 'HPAI', or 'check seq' depending on the number of basic residues.

    > **NOTE:** Only sequences with sequence headers ending with 'HA1' or 'HA2' are processed, e.g.
    `>SAMPLE|other_info HA1`. If either HA1 or HA2 is missing, the sequence is skipped.

    ## Output

    The results are written to a TSV file with the following columns:

    * **Sample**: Sample ID

    * **Cleavage Sequence**: Detected cleavage sequence (HA1 motif matching sequence + first 3 HA2 residues)

    * **Residue Type**: Classification of the basic residues in the HA1 motif

    * **Basic Residue Count**: Number of basic residues in the HA1 motif

    * **Pathogenicity**: Classification of the pathogenicity based on the basic residue count

    * HA1 Subsequence (Last 15 residues)

    * HA2 Subsequence (First 3 residues)

    * HA1 Sequence Header

    * HA1 Amino Acid Sequence

    * HA2 Sequence Header

    * HA2 Amino Acid Sequence
    """
    init_logging(verbose)
    if output_tsv.exists() and not force:
        logger.error(f"Output file '{output_tsv}' exists. Use --force to overwrite.")
        raise typer.Abort()
    logger.info(f"{input_fasta=}")
    # Step 1: Extract HA1 and HA2 sequences
    ha1_ha2_seqs = get_ha1_ha2_seqs(input_fasta)
    logger.info(f"Read {len(ha1_ha2_seqs)} HA1/HA2 sequences from '{input_fasta}'")
    logger.debug(f"{ha1_ha2_seqs=}")
    # Step 2: Process cleavage site analysis and classification
    sites = find_cleavage_sites(ha1_ha2_seqs)
    write_sites(output_tsv, sites)


if __name__ == "__main__":
    app()
