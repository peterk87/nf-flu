#!/usr/bin/env python

import logging
from pathlib import Path

import typer
from rich.logging import RichHandler


def init_logging():
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def main(
        input_fasta: Path = typer.Option(..., "--input-fasta", "-f"),
        input_feature_table: Path = typer.Option(..., "--input-feature-table", "-t"),
        outdir: Path = typer.Option(Path("./"), "--outdir", "-o"),
        prefix: str = typer.Option("SAMPLE", "--prefix", "-p"),
):
    init_logging()
    logging.info(f"{input_fasta=}")
    logging.info(f"{input_feature_table=}")
    logging.info(f"{outdir=}")
    logging.info(f"{prefix=}")
    ft_seqids = []
    ft_out = outdir / f"{prefix}.tbl"
    with input_feature_table.open() as fh, ft_out.open("w") as fout:
        seqid_count = 0
        for line in fh:
            if line.startswith(">Feature "):
                seqid = line.replace(">Feature ", "").strip()
                seqid_count += 1
                new_seqid = f"SEQUENCE-{seqid_count}"
                ft_seqids.append((seqid, new_seqid))
            fout.write(line.replace(seqid, new_seqid))
    logging.info(f"table2asn compatible feature table written to '{ft_out}'")

    namesub_path = outdir / f"{prefix}.namesub.txt"
    with namesub_path.open("w") as fout:
        for seqid, new_seqid in ft_seqids:
            fout.write(f"{new_seqid}\t{seqid}\n")
    logging.info(
        f"Text file with original and table2asn seqids written to '{namesub_path}'",
    )
    fasta_out = outdir / f"{prefix}.fa"
    seqid_to_new_seqid = dict(ft_seqids)
    logging.info(f"{seqid_to_new_seqid=}")
    with input_fasta.open() as fh, fasta_out.open("w") as fout:
        for line in fh:
            if line.startswith(">"):
                seqid, *_ = line[1:].strip().split(" ")
                try:
                    fout.write(line.replace(seqid, seqid_to_new_seqid[seqid]))
                except KeyError:
                    logging.warning(f"{seqid=} not in feature table seqids!")
            else:
                fout.write(line)
    logging.info(f"table2asn compatible fasta written to '{fasta_out}'")


if __name__ == "__main__":
    typer.run(main)
