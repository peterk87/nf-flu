#!/usr/bin/env python

"""
Generate an Influenza H/N subtyping report from nucleotide BLAST results for one or more genomes.
Metadata from the NCBI Influenza DB is merged in with the BLAST results to improve subtyping results and provide context for the results obtained.
For reference:
Segment 4 - hemagglutinin (HA) gene
Segment 6 - neuraminidase (NA) gene
"""

import logging
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import click
import numpy as np
import pandas as pd
import polars as pl
from rich.logging import RichHandler

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

pl.enable_string_cache(True)

influenza_segment = {
    1: "1_PB2",
    2: "2_PB1",
    3: "3_PA",
    4: "4_HA",
    5: "5_NP",
    6: "6_NA",
    7: "7_M",
    8: "8_NS",
}

# Column names/types/final report names
blast_cols = [
    ("qaccver", str),
    ("saccver", str),
    ("pident", float),
    ("length", pl.UInt16),
    ("mismatch", pl.UInt16),
    ("gapopen", pl.UInt16),
    ("qstart", pl.UInt16),
    ("qend", pl.UInt16),
    ("sstart", pl.UInt16),
    ("send", pl.UInt16),
    ("evalue", pl.Float32),
    ("bitscore", pl.Float32),
    ("qlen", pl.UInt16),
    ("slen", pl.UInt16),
    ("qcovs", pl.Float32),
    ("stitle", str),
]

blast_results_report_columns = [
    ("sample", "Sample"),
    ("sample_segment", "Sample Genome Segment Number"),
    ("#Accession", "Reference NCBI Accession"),
    ("Genotype", "Reference Subtype"),
    ("pident", "BLASTN Percent Identity"),
    ("length", "BLASTN Alignment Length"),
    ("mismatch", "BLASTN Mismatches"),
    ("gapopen", "BLASTN Gaps"),
    ("qstart", "BLASTN Sample Start Index"),
    ("qend", "BLASTN Sample End Index"),
    ("sstart", "BLASTN Reference Start Index"),
    ("send", "BLASTN Reference End Index"),
    ("evalue", "BLASTN E-value"),
    ("bitscore", "BLASTN Bitscore"),
    ("qlen", "Sample Sequence Length"),
    ("slen", "Reference Sequence Length"),
    ("qcovs", "Sample Sequence Coverage of Reference Sequence"),
    ("stitle", "Reference Sequence ID"),
    ("Segment", "Reference Genome Segment Number"),
    ("GenBank_Title", "Reference Virus Name"),
    ("Host", "Reference Host"),
    ("Geo_Location", "Reference Geo Location"),
    ("Collection_Date", "Reference Collection Date"),
    ("Release_Date", "Reference Release Date"),
]

subtype_results_summary_columns = [
    "sample",
    "Genotype",
    "H_top_accession",
    "H_type",
    "H_virus_name",
    "H_NCBI_Influenza_DB_proportion_matches",
    "N_top_accession",
    "N_type",
    "N_virus_name",
    "N_NCBI_Influenza_DB_proportion_matches",
]

columns_H_summary_results = [
    "sample",
    "Genotype",
    "H_top_accession",
    "H_NCBI_Influenza_DB_proportion_matches",
    "H_NCBI_Influenza_DB_subtype_matches",
    "H_NCBI_Influenza_DB_total_matches",
    "H_sample_segment_length",
    "H_top_align_length",
    "H_top_bitscore",
    "H_top_country",
    "H_top_date",
    "H_top_gaps",
    "H_top_host",
    "H_top_mismatch",
    "H_top_pident",
    "H_top_seq_length",
    "H_type",
    "H_virus_name",
]

columns_N_summary_results = [
    "sample",
    "Genotype",
    "N_top_accession",
    "N_NCBI_Influenza_DB_proportion_matches",
    "N_NCBI_Influenza_DB_subtype_matches",
    "N_NCBI_Influenza_DB_total_matches",
    "N_sample_segment_length",
    "N_top_align_length",
    "N_top_bitscore",
    "N_top_country",
    "N_top_date",
    "N_top_gaps",
    "N_top_host",
    "N_top_mismatch",
    "N_top_pident",
    "N_top_seq_length",
    "N_type",
    "N_virus_name",
]

subtype_results_summary_final_names = {
    "sample": "Sample",
    "Genotype": "Subtype Prediction",
    "N_type": "N: type prediction",
    "N_top_accession": "N: top match accession",
    "N_virus_name": "N: top match virus name",
    "N_top_host": "N: top match host",
    "N_top_date": "N: top match collection date",
    "N_top_country": "N: top match country",
    "N_top_pident": "N: top match BLASTN % identity",
    "N_top_align_length": "N: top match BLASTN alignment length",
    "N_top_mismatch": "N: top match BLASTN mismatches",
    "N_top_gaps": "N: top match BLASTN gaps",
    "N_top_bitscore": "N: top match BLASTN bitscore",
    "N_top_seq_length": "N: top match sequence length",
    "N_sample_segment_length": "N: sample segment length",
    "N_NCBI_Influenza_DB_proportion_matches": "N: NCBI Influenza DB subtype match proportion",
    "N_NCBI_Influenza_DB_subtype_matches": "N: NCBI Influenza DB subtype match count",
    "N_NCBI_Influenza_DB_total_matches": "N: NCBI Influenza DB total count",
    "H_type": "H: type prediction",
    "H_top_accession": "H: top match accession",
    "H_virus_name": "H: top match virus name",
    "H_top_host": "H: top match host",
    "H_top_date": "H: top match collection date",
    "H_top_country": "H: top match country",
    "H_top_pident": "H: top match BLASTN % identity",
    "H_top_align_length": "H: top match BLASTN alignment length",
    "H_top_mismatch": "H: top match BLASTN mismatches",
    "H_top_gaps": "H: top match BLASTN gaps",
    "H_top_bitscore": "H: top match BLASTN bitscore",
    "H_top_seq_length": "H: top match sequence length",
    "H_sample_segment_length": "H: sample segment length",
    "H_NCBI_Influenza_DB_proportion_matches": "H: NCBI Influenza DB subtype match proportion",
    "H_NCBI_Influenza_DB_subtype_matches": "H: NCBI Influenza DB subtype match count",
    "H_NCBI_Influenza_DB_total_matches": "H: NCBI Influenza DB total count",
}

# Regex to find unallowed characters in Excel worksheet names
REGEX_UNALLOWED_EXCEL_WS_CHARS = re.compile(r"[\\:/?*\[\]]+")


def parse_blast_result(
        blast_result: str,
        df_metadata: pl.DataFrame,
        regex_subtype_pattern: str,
        get_top_ref: bool,
        top: int = 3,
        pident_threshold: float = 0.85,
        min_aln_length: int = 50,
) -> Optional[Tuple[pl.DataFrame, Dict]]:
    logging.info(f"Parsing BLAST results from {blast_result}")

    try:
        df_filtered = (
            pl.scan_csv(
                blast_result,
                has_header=False,
                separator="\t",
                new_columns=[name for name, coltype in blast_cols],
                dtypes=dict(blast_cols),
            )
            .filter(
                (pl.col("pident") >= (pident_threshold * 100))
                & (pl.col("length") >= min_aln_length)
            )
            .collect(streaming=True)
        )
    except pl.exceptions.NoDataError:
        logging.warning(f"No BLAST results found in {blast_result}")
        return None

    sample_name: str = re.sub(r"^(.+)_\d$", r"\1", df_filtered["qaccver"][0])
    logging.info(
        f"{sample_name} | n={df_filtered.shape[0]} | Filtered for hits above {pident_threshold}% identity."
        f"and Min Alignment length > {min_aln_length}"
    )
    df_filtered = df_filtered.with_columns([
        pl.col('saccver').str.strip().alias("#Accession"),
        pl.lit(sample_name, dtype=pl.Categorical).alias("sample"),
        pl.col('qaccver').str.extract(r".+_(\d)$").cast(pl.Categorical).alias("sample_segment"),
        pl.col("stitle").str.extract(regex_subtype_pattern).alias("subtype_from_match_title").cast(pl.Categorical)
    ])
    logging.info(
        f"{sample_name} | Merging NCBI Influenza DB genome metadata with BLAST results on accession."
    )
    df_merge = df_filtered.join(df_metadata, on="#Accession", how="left")
    del df_filtered
    del df_metadata
    df_merge = df_merge.with_columns(
        pl.when(pl.col("Genotype").is_null())
        .then(pl.col("subtype_from_match_title"))
        .otherwise(pl.col("Genotype"))
        .alias("Genotype")
    )
    df_merge = df_merge.sort(
        by=["sample_segment", "bitscore"], descending=[False, True]
    )

    segments = df_merge["sample_segment"].unique().sort()
    dfs = [
        df_merge.filter(pl.col("sample_segment") == seg).head(top)
        for seg in segments
    ]
    df_top_seg_matches = pl.concat(dfs, how="vertical")
    cols = pl.Series([x for x, _ in blast_results_report_columns])
    df_top_seg_matches = df_top_seg_matches.select(pl.col(cols))
    subtype_results_summary = {"sample": sample_name}
    if not get_top_ref:
        is_iav = not df_top_seg_matches.select(pl.col("Genotype").is_null().all())[0, 0]
        H_results = None
        N_results = None
        if "4" in segments:
            H_results = find_h_or_n_type(df_merge, "4", is_iav)
            subtype_results_summary |= H_results
        if "6" in segments:
            N_results = find_h_or_n_type(df_merge, "6", is_iav)
            subtype_results_summary.update(N_results)
        subtype_results_summary["Genotype"] = get_subtype_value(H_results, N_results, is_iav)

    return df_top_seg_matches, subtype_results_summary


def get_subtype_value(H_results: Optional[Dict], N_results: Optional[Dict], is_iav: bool) -> str:
    subtype = ""
    if not is_iav:
        return "N/A"
    if H_results is None and N_results is None:
        subtype = "-"
    elif H_results is not None and N_results is None:
        H: str = H_results.get("H_type", "")
        subtype = f"H{H}" if H != "" else "-"
    elif H_results is None:
        N: str = N_results.get("N_type", "")
        subtype = f"N{N}" if N != "" else "-"
    else:
        H: str = H_results.get("H_type", "")
        N: str = N_results.get("N_type", "")
        if H or N:
            if H != "":
                H = f"H{H}"
            if N != "":
                N = f"N{N}"
            subtype = f"{H}{N}"
        else:
            subtype = "-"
    return subtype


def find_h_or_n_type(df_merge, seg, is_iav):
    assert seg in [
        "4",
        "6",
    ], "Can only determine H or N type from segments 4 or 6, respectively!"
    type_name = "H_type" if seg == "4" else "N_type"
    h_or_n = type_name[0]
    df_segment = df_merge.filter(pl.col("sample_segment") == seg)
    if is_iav:
        type_counts = df_segment["Genotype"].value_counts(sort=True)
        type_counts = type_counts.filter(~pl.col("Genotype").is_null())
        df_type_counts = type_counts.with_columns(pl.lit(type_counts["Genotype"].str.extract(reg_h_or_n_type + r"(\d+)").alias(type_name)))
        df_type_counts = df_type_counts.filter(~pl.col(type_name).is_null())
        logging.debug(f"{df_type_counts}")
        type_to_count = defaultdict(int)
        for x in df_type_counts.iter_rows(named=True):
            type_to_count[x[type_name]] += x["counts"]
        type_to_count = list(type_to_count.items())
        type_to_count.sort(key=lambda x: x[1], reverse=True)
        top_type, top_type_count = type_to_count[0]
        total_count = type_counts["counts"].sum()
        logging.info(
            f"{h_or_n}{top_type} n={top_type_count}/{total_count} ({top_type_count / total_count:.1%})"
        )
        df_segment = df_segment.with_columns(
            pl.lit(
                df_segment["Genotype"]
                .str.contains(f".*{reg_h_or_n_type}" + top_type + r".*")
                .fill_null(False)
                .alias("type_mask")
            )
        )
        df_seg_top_type = df_segment.filter(pl.col("type_mask") == True).drop("type_mask")
        top_result: pl.Series = list(df_seg_top_type.head(1).iter_rows(named=True))[0]
    else:
        top_type = "N/A"
        top_type_count = "N/A"
        total_count = "N/A"
        top_result: pl.Series = list(df_segment.head(1).iter_rows(named=True))[0]

    results_summary = {
        f"{h_or_n}_type": top_type if is_iav else "N/A",
        f"{h_or_n}_sample_segment_length": top_result["qlen"],
        f"{h_or_n}_top_pident": top_result["pident"],
        f"{h_or_n}_top_mismatch": top_result["mismatch"],
        f"{h_or_n}_top_gaps": top_result["gapopen"],
        f"{h_or_n}_top_bitscore": top_result["bitscore"],
        f"{h_or_n}_top_align_length": top_result["length"],
        f"{h_or_n}_top_accession": top_result["#Accession"],
        f"{h_or_n}_top_host": top_result["Host"],
        f"{h_or_n}_top_country": top_result["Geo_Location"],
        f"{h_or_n}_top_date": top_result["Collection_Date"],
        f"{h_or_n}_top_seq_length": top_result["slen"],
        f"{h_or_n}_virus_name": top_result["GenBank_Title"],
        f"{h_or_n}_NCBI_Influenza_DB_subtype_matches": top_type_count,
        f"{h_or_n}_NCBI_Influenza_DB_total_matches": total_count,
        f"{h_or_n}_NCBI_Influenza_DB_proportion_matches": top_type_count / total_count if is_iav else "N/A",
    }
    logging.info(f"Seg {seg} results: {results_summary}")
    return results_summary


@click.command()
@click.option("-m", "--flu-metadata", help="NCBI Influenza genomeset.dat metadata file")
@click.option("-x", "--excel-report", default="report.xlsx", help="Excel report")
@click.option("--top", default=3, help="Top N matches to each segment to report")
@click.option(
    "--pident-threshold", default=0.85, help="BLAST percent identity threshold"
)
@click.option('--min-aln-length', default=50, help="Min BLAST alignment length threshold")
@click.option("--threads", default=4, help="Number of BLAST result parsing threads.")
@click.option("--get-top-ref", default=False, help="Get top ref accession id from ncbi database.")
@click.option("--sample-name", default="", help="Sample Name.")
@click.argument("blast_results", nargs=-1)
def report(flu_metadata, blast_results, excel_report, top, pident_threshold,
           min_aln_length, threads, get_top_ref, sample_name):
    from rich.traceback import install
    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    logging.info(f'Parsing Influenza metadata file "{flu_metadata}"')

    md_cols = [
        ("#Accession", str),
        ("Release_Date", pl.Categorical),
        ("Genus", pl.Categorical),
        ("Length", pl.UInt16),
        ("Genotype", str),
        ("Segment", pl.Categorical),
        ("Publications", str),
        ("Geo_Location", pl.Categorical),
        ("Host", pl.Categorical),
        ("Isolation_Source", pl.Categorical),
        ("Collection_Date", pl.Categorical),
        ("GenBank_Title", str),
    ]
    df_md = pl.read_csv(
        flu_metadata,
        has_header=True,
        has_header=False,
        dtypes=dict(md_cols),
    )

    unique_subtypes = df_md.select("Genotype").unique()
    unique_subtypes = unique_subtypes.filter(~pl.col("Genotype").is_null())
    logging.info(
        f"Parsed Influenza metadata file into DataFrame with n={df_md.shape[0]} rows and n={df_md.shape[1]} columns. There are {len(unique_subtypes)} unique subtypes. "
    )
    regex_subtype_pattern = r"\((H\d+N\d+|" + "|".join(list(unique_subtypes["Genotype"])) + r")\)"
    results = [
        parse_blast_result(blast_result, df_md, regex_subtype_pattern, get_top_ref, top=top,
                           pident_threshold=pident_threshold,
                           min_aln_length=min_aln_length) for blast_result in blast_results]

    if not get_top_ref:
        dfs_blast = []
        all_subtype_results = {}
        for parsed_result in results:
            if parsed_result is None:
                continue
            df_blast, subtype_results_summary = parsed_result
            if df_blast is not None:
                dfs_blast.append(df_blast.to_pandas())
            sample = subtype_results_summary["sample"]
            all_subtype_results[sample] = subtype_results_summary
        df_all_blast = pd.concat(dfs_blast).rename(
            columns=dict(blast_results_report_columns)
        )
        df_subtype_results = pd.DataFrame(all_subtype_results).transpose()
        cols = pd.Series(subtype_results_summary_columns)
        cols = cols[cols.isin(df_subtype_results.columns)]
        df_subtype_predictions = df_subtype_results[cols].rename(
            columns=subtype_results_summary_final_names
        )
        cols = pd.Series(columns_H_summary_results)
        cols = cols[cols.isin(df_subtype_results.columns)]
        df_H = df_subtype_results[cols].rename(columns=subtype_results_summary_final_names)
        cols = pd.Series(columns_N_summary_results)
        cols = cols[cols.isin(df_subtype_results.columns)]
        df_N = df_subtype_results[cols].rename(columns=subtype_results_summary_final_names)
        # Add segment name for more informative
        df_all_blast["Sample Genome Segment Number"] = df_all_blast["Sample Genome Segment Number"]. \
            apply(lambda x: influenza_segment[int(x)])
        write_excel(
            [
                ("Subtype Predictions", df_subtype_predictions),
                ("Top Segment Matches", df_all_blast),
                ("H Segment Results", df_H),
                ("N Segment Results", df_N),
            ],
            output_dest=excel_report,
        )
    else:
        df_blast, subtype_results_summary = results[0]
        df_blast = df_blast.rename(mapping=dict(blast_results_report_columns))
        df_ref_id = df_blast.select(
            pl.col([
                'Sample',
                'Sample Genome Segment Number',
                'Reference NCBI Accession',
                'BLASTN Bitscore',
                'Reference Sequence ID'
            ])
        )
        df_ref_id = df_ref_id.with_columns(
            pl.when(pl.col("Reference NCBI Accession").is_null())
            .then(pl.col("Reference Sequence ID"))
            .otherwise(pl.col("Reference NCBI Accession"))
            .str.strip()
            .alias('Reference NCBI Accession')
        )
        df_ref_id = df_ref_id.with_columns(
            pl.col("Sample Genome Segment Number").apply(lambda x: influenza_segment[int(x)])
            .alias("Sample Genome Segment Number"))
        df_ref_id.write_csv(sample_name + ".topsegments.csv", separator=",", has_header=True)


def get_col_widths(df, index=False):
    """Calculate column widths based on column headers and contents"""
    if index:
        yield max(
            [len(str(s)) for s in df.index.values] + [len(str(df.index.name))]
        )
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
                fixed_sheetname = f"{idx}_{fixed_sheetname[:max_chars]}"
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

            logging.info(f'Writing table to Excel sheet "{fixed_sheetname}"')
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
