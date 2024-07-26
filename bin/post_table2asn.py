#!/usr/bin/env python

import logging
import os
import re
from pathlib import Path

import typer
from Bio import SeqIO
from Bio.SeqIO.InsdcIO import _insdc_location_string
from Bio.SeqRecord import SeqRecord
from gfflu.io import get_translation, write_gff
from rich.logging import RichHandler

REGEX_SUBNAME = re.compile(r"(SEQUENCE-\d+)")


def init_logging():
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def get_namesub(namesub_txt: Path) -> dict[str, str]:
    namesub = {}
    with open(namesub_txt) as fh:
        for line in fh:
            line = line.strip()
            if line == "":
                continue
            sub_name, original_name = line.split("\t", maxsplit=1)
            namesub[sub_name] = original_name
    return namesub


def write_aa_fasta(recs: list[SeqRecord], out_file: os.PathLike) -> None:
    """Write out FASTA file with amino acid sequences"""
    with open(out_file, "w") as out_handle:
        for rec in recs:
            for feature in rec.features:
                if feature.type not in {"CDS", "mat_peptide", "sig_peptide"}:
                    continue
                try:
                    translation = feature.qualifiers["translation"][0]
                except (KeyError, IndexError):
                    logging.warning(
                        f"Could not retrieve translation qualifier for feature {feature} from {rec}. "
                        "Trying to obtain translation from nucleotide sequence."
                    )
                    translation = get_translation(rec.seq, feature.location)
                try:
                    gene = feature.qualifiers["gene"][0]
                except KeyError:
                    gene = ""
                try:
                    product = feature.qualifiers["product"][0]
                except KeyError:
                    product = ""
                location_str = _insdc_location_string(
                    location=feature.location, rec_length=len(rec.seq)
                )
                out_handle.write(
                    f">{rec.id}|{feature.type}|{gene}|{location_str} {product}\n{translation}\n"
                )


def write_cds_nt_fasta(recs: list[SeqRecord], out_file: os.PathLike) -> None:
    """Write out FASTA file with CDS nt sequences"""
    with open(out_file, "w") as out_handle:
        for rec in recs:
            for feature in rec.features:
                if feature.type not in {"CDS", "mat_peptide", "sig_peptide"}:
                    continue
                subseq = str(feature.location.extract(rec.seq))
                try:
                    gene = feature.qualifiers["gene"][0]
                except KeyError:
                    gene = ""
                try:
                    product = feature.qualifiers["product"][0]
                except KeyError:
                    product = ""
                location_str = _insdc_location_string(
                    location=feature.location, rec_length=len(rec.seq)
                )
                out_handle.write(
                    f">{rec.id}|{feature.type}|{gene}|{location_str} {product}\n{subseq}\n"
                )


def output_subbed_genbank(
    gbk_path: Path,
    gbk_out: Path,
    namesub: dict[str, str],
) -> None:
    with open(gbk_path) as fh, open(gbk_out, "w") as fout:
        for line in fh:
            match = REGEX_SUBNAME.search(line)
            if match:
                subname = match.group(1)
                original_name = namesub[subname]
                line = line.replace(subname, original_name)
            fout.write(line)


def main(
    table2asn_genbank: Path,
    namesub_txt: Path,
    outdir: Path = typer.Option(Path("./"), "--outdir", "-o"),
    prefix: str = typer.Option("SAMPLE", "--prefix", "-p"),
):
    init_logging()
    logging.info(f"{table2asn_genbank=}")
    logging.info(f"{namesub_txt=}")
    logging.info(f"{outdir=}")
    logging.info(f"{prefix=}")
    namesub = get_namesub(namesub_txt)
    logging.info(f'Read {len(namesub)} from "{namesub_txt}"')
    gbk_out = outdir / f"{prefix}.gbk"
    output_subbed_genbank(table2asn_genbank, gbk_out, namesub)
    logging.info(f'Wrote Genbank with original names to "{gbk_out}"')
    recs = list(SeqIO.parse(gbk_out, "genbank"))
    logging.info(
        f"Read {len(recs)} sequence records from '{gbk_out}'. Writing GFF, amino acid and nucleotide FASTAs for "
        "sequence features (CDS, mature peptides, signal peptides)..."
    )
    gff_out = outdir / f"{prefix}.gff"
    write_gff(recs, gff_out)
    logging.info(f'Wrote GFF to "{gff_out}"')
    faa_out = outdir / f"{prefix}.faa"
    write_aa_fasta(recs, faa_out)
    logging.info(f'Wrote amino acid FASTA to "{faa_out}"')
    ffn_out = outdir / f"{prefix}.ffn"
    write_cds_nt_fasta(recs, ffn_out)
    logging.info(f'Wrote nucleotide FASTA of gene features to "{faa_out}"')
    logging.info("Done!")


if __name__ == "__main__":
    typer.run(main)
