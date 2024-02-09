#!/usr/bin/env python

from typing import List, Tuple
import re
import logging

import click
import pandas as pd
import numpy as np

from rich.logging import RichHandler

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

# Column names/types/final report names
blast_cols = [
    ("qaccver", "category"),
    ("saccver", str),
    ("pident", float),
    ("length", "uint16"),
    ("mismatch", "uint16"),
    ("gapopen", "uint16"),
    ("qstart", "uint16"),
    ("qend", "uint16"),
    ("sstart", "uint16"),
    ("send", "uint16"),
    ("evalue", np.float16),
    ("bitscore", np.float16),
    ("qlen", "uint16"),
    ("slen", "uint16"),
    ("qcovs", np.float16),
    ("stitle", str),
]

# Regex to find unallowed characters in Excel worksheet names
REGEX_UNALLOWED_EXCEL_WS_CHARS = re.compile(r"[\\:/?*\[\]]+")


@click.command()
@click.option("-x", "--excel-report", default="report.xlsx", help="Excel report")
@click.option('--min-aln-length', default=50, help="Min BLAST alignment length threshold")
@click.option("-b", "--blast_results", default="", help="Blast Result.")
def report(blast_results, excel_report, min_aln_length):
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    df_blast_result = pd.read_csv(
        blast_results,
        sep="\t",
        names=[name for name, coltype in blast_cols],
        dtype={name: coltype for name, coltype in blast_cols},
    )
    if df_blast_result.empty:
        # write empty report, to keep workflow run without stopping due this error
        write_excel(
            [
                ("Mismatch_Report", df_blast_result),
                ("Blastn_Results", df_blast_result)
            ],
            output_dest=excel_report,
        )
        logging.warning(f"No BLASTN results in {blast_results}!, empty report")
        return
    df_blast_result["ref_name"] = df_blast_result["stitle"].str.extract('(.+?)_[sS]egment')
    df_blast_result["sample_name"] = df_blast_result["qaccver"].str.extract('(.+?)_[1-9]$')
    df_blast_result["segment_name"] = df_blast_result["qaccver"].str.extract(r".+_(\d)$")
    df_blast_result.drop(columns=["qaccver", "saccver", "stitle"], inplace=True)
    df_blast_result = df_blast_result.reindex(
        columns=["sample_name", "segment_name", "ref_name", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"])

    # get df for mismatch report
    df_filtered = df_blast_result[df_blast_result["length"] >= min_aln_length]
    segments = df_blast_result["segment_name"].unique()
    ref_names = df_blast_result["ref_name"].unique()
    df_mismatch_report = pd.DataFrame(index=segments, columns=ref_names)
    for segment in segments:
        for ref_name in ref_names:
            mismatch = df_filtered.query("ref_name == @ref_name and segment_name == @segment")["mismatch"].values
            if len(mismatch):
                df_mismatch_report.loc[segment, ref_name] = mismatch[0]
            else:
                df_mismatch_report.loc[segment, ref_name] = ''
    df_mismatch_report.insert(0, "Segment", segments)
    df_mismatch_report.loc["Total"] = pd.Series(df_mismatch_report[ref_names].sum())
    df_blast_result.columns = ["Sample", "Sample Genome Segment Number", "Reference Virus Name",
                               "BLASTN Percent Identity",
                               "BLASTN Alignment Length", "BLASTN Mismatches", "BLASTN Gaps",
                               "BLASTN Sample Start Index",
                               "BLASTN Sample End Index", "BLASTN Reference Start Index", "BLASTN Reference End Index",
                               "BLASTN E-value", "BLASTN Bitscore", "Sample Sequence Length",
                               "Reference Sequence Length",
                               "Sample Sequence Coverage of Reference Sequence"]
    write_excel(
        [
            ("Mismatch_Report", df_mismatch_report),
            ("Blastn_Results", df_blast_result)
        ],
        output_dest=excel_report,
    )


def get_col_widths(df, index=False):
    """Calculate column widths based on column headers and contents"""
    if index:
        idx_max = max(
            [len(str(s)) for s in df.index.values] + [len(str(df.index.name))]
        )
        yield idx_max
    for c in df.columns:
        # get max length of column contents and length of column header
        yield np.max([df[c].astype(str).str.len().max() + 1, len(c) + 1])


def write_excel(
        name_dfs: List[Tuple[str, pd.DataFrame]],
        output_dest: str,
        sheet_name_index: bool = True,
) -> None:
    logging.info("Starting to write tabular data to worksheets in Excel workbook")
    with pd.ExcelWriter(output_dest, engine="xlsxwriter") as writer:
        idx = 1
        for name_df in name_dfs:
            if not isinstance(name_df, (list, tuple)):
                logging.error(
                    'Input "%s" is not a list or tuple (type="%s"). Skipping...',
                    name_df,
                    type(name_df),
                )
                continue
            sheetname, df = name_df
            fixed_sheetname = REGEX_UNALLOWED_EXCEL_WS_CHARS.sub("_", sheetname)
            # fixed max number of characters in sheet name due to compatibility
            if sheet_name_index:
                max_chars = 28
                fixed_sheetname = "{}_{}".format(idx, fixed_sheetname[:max_chars])
            else:
                max_chars = 31
                fixed_sheetname = fixed_sheetname[:max_chars]

            if len(fixed_sheetname) > max_chars:
                logging.warning(
                    'Sheetname "%s" is >= %s characters so may be truncated (n=%s)',
                    max_chars,
                    fixed_sheetname,
                    len(fixed_sheetname),
                )

            logging.info('Writing table to Excel sheet "{}"'.format(fixed_sheetname))
            df.to_excel(
                writer, sheet_name=fixed_sheetname, index=False, freeze_panes=(1, 1)
            )
            worksheet = writer.book.get_worksheet_by_name(fixed_sheetname)
            for i, width in enumerate(get_col_widths(df, index=False)):
                worksheet.set_column(i, i, width)
            idx += 1
    logging.info('Done writing worksheets to spreadsheet "%s".', output_dest)


if __name__ == "__main__":
    report()
