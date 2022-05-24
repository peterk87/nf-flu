#!/usr/bin/env python

from pathlib import Path

import pandas as pd
import typer
from rich.logging import RichHandler
import logging


def check_sample_names(df: pd.DataFrame) -> None:
    have_whitespace = df['sample'].str.contains(r'\s', regex=True)
    n_samples_with_whitespace = have_whitespace.sum()
    if have_whitespace.sum() > 0:
        raise ValueError(
            f'Found {n_samples_with_whitespace} sample names with whitespace: '
            f'{"; ".join(df.loc[have_whitespace, "sample"])}\n'
            f'{df.loc[have_whitespace, :]}\n'
            f'Please check your sample sheet. Sample names should not have spaces.'
        )


def adjust_reads_path(p: str) -> str:
    assert (
            p.endswith(".fastq")
            or p.endswith(".fastq.gz")
            or p.endswith(".fq")
            or p.endswith(".fq.gz")
    ), 'FASTQ file "{p}" does not have expected extension: ".fastq", ".fastq.gz", ".fq", ".fq.gz"'
    if p.startswith("http") or p.startswith("ftp"):
        return p
    else:
        path = Path(p)
        return str(path.resolve().absolute())


def main(input_path: Path, platform: str, output_sample_sheet: Path):
    """Check and reformat sample sheet into CSV

    Outputs CSV with headers: sample, fastq1, fastq2, single_end"""
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    logging.info(
        f'input_path="{input_path}" output_sample_sheet="{output_sample_sheet}"'
    )
    ext = input_path.suffix.lower()
    logging.info(f"Input sample sheet extension: {ext}")
    try:
        if ext in [".tsv", ".txt", ".tab"]:
            df = pd.read_table(input_path, dtype="str")
        elif ext == ".csv":
            df = pd.read_csv(input_path, dtype="str")
        elif ext in [".xls", ".xlsx", ".ods"]:
            df = pd.read_excel(input_path, dtype="str")
        else:
            raise ValueError(f'Unknown file format for sample sheet "{input_path}"')
    except Exception as ex:
        logging.exception(ex)
        raise ex
    if platform == 'illumina':
        assert (
                df.shape[1] == 3
        ), f"3 columns expected in sample sheet, but {df.shape[1]} found!"
        df.columns = ["sample", "fastq1", "fastq2"]
        check_sample_names(df)
        fastq1_paths = []
        fastq2_paths = []
        single_ends = []
        for row in df.itertuples():
            fastq1 = row.fastq1
            fastq2 = row.fastq2
            fastq1_isnull = pd.isnull(fastq1)
            fastq2_isnull = pd.isnull(fastq2)
            if not fastq1_isnull and not fastq2_isnull:
                single_ends.append(False)
                fastq1_paths.append(adjust_reads_path(fastq1))
                fastq2_paths.append(adjust_reads_path(fastq2))
            elif not fastq1_isnull and fastq2_isnull:
                single_ends.append(True)
                fastq1_paths.append(adjust_reads_path(fastq1))
                fastq2_paths.append(None)
            elif fastq1_isnull and not fastq2_isnull:
                single_ends.append(True)
                fastq1_paths.append(adjust_reads_path(fastq2))
                fastq2_paths.append(None)
            else:
                err_msg = f'Forward and/or reverse reads paths NOT specified for sample "{row.sample}"'
                raise ValueError(err_msg)

        df["fastq1"] = fastq1_paths
        df["fastq2"] = fastq2_paths
        df["single_end"] = single_ends
        df.to_csv(output_sample_sheet, index=False)
        logging.info(f'Wrote reformatted sample sheet CSV to "{output_sample_sheet}"')
    elif platform == 'nanopore':
        assert (
                df.shape[1] == 2
        ), f"2 columns expected in sample sheet, but {df.shape[1]} found!"
        df.columns = ["sample", "barcode"]
        check_sample_names(df)
        df.to_csv(output_sample_sheet, index=False)
        logging.info(f'Wrote reformatted sample sheet CSV to "{output_sample_sheet}"')


if __name__ == "__main__":
    typer.run(main)
