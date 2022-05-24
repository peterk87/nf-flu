#!/usr/bin/env python

"""
Generate an Influenza H/N subtyping report from nucleotide BLAST results for one or more genomes.
Metadata from the NCBI Influenza DB is merged in with the BLAST results to improve subtyping results and provide context for the results obtained.
For reference:
Segment 4 - hemagglutinin (HA) gene
Segment 6 - neuraminidase (NA) gene
"""

from typing import Dict, List, Optional, Tuple
import re
import logging
from collections import defaultdict
from multiprocessing import Pool

import click
import pandas as pd
import numpy as np
from rich.console import Console
from rich.logging import RichHandler

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

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

blast_results_report_columns = [
    ("sample", "Sample"),
    ("sample_segment", "Sample Genome Segment Number"),
    ("accession", "Reference NCBI Accession"),
    ("subtype", "Reference Subtype"),
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
    ("segment", "Reference Genome Segment Number"),
    ("virus_name", "Reference Virus Name"),
    ("host", "Reference Host"),
    ("country", "Reference Country"),
    ("date", "Reference Collection Date"),
    ("age", "Reference Patient Age"),
    ("gender", "Reference Patient Gender"),
    ("group_id", "Reference Group ID"),
]

subtype_results_summary_columns = [
    "sample",
    "subtype",
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
    "subtype",
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
    "subtype",
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
    "subtype": "Subtype Prediction",
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
        df_metadata: pd.DataFrame,
        regex_subtype_pattern: str,
        top: int = 3,
        pident_threshold: float = 0.85,
        min_aln_length: int = 50,
) -> Optional[
    Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame], Dict]
]:
    logging.info(f"Parsing BLAST results from {blast_result}")

    df = pd.read_csv(
        blast_result,
        sep="\t",
        names=[name for name, coltype in blast_cols],
        dtype={name: coltype for name, coltype in blast_cols},
    )
    if df.empty:
        logging.error(f"No BLASTN results in {blast_result}!")
        return
    # Assuming that all BLAST result entries have the same sample name so
    # extracting the sample name from the first entry qaccver field
    sample_name: str = re.sub(r"^(.+)_\d$", r"\1", df["qaccver"][0])
    logging.info(f"Parsed {df.shape[0]} BLAST results from {blast_result}")
    logging.info(
        f"{sample_name} | n={df.shape[0]} | Filtering for hits above {pident_threshold}% identity."
    )
    df_filtered = df[
        (df["pident"] >= (pident_threshold * 100)) & (df["length"] >= min_aln_length)
        ]
    logging.info(
        f"{sample_name} | n={df_filtered.shape[0]} | Filtered for hits above {pident_threshold}% identity."
    )
    # Sequences header has been corrected by GUNZIP Process
    # zcat $archive | sed -E 's/^>gi\\|[0-9]+\\|gb\\|(\\w+)\\|(.*)/>\\1 \\2/' > influenza.fna
    df_filtered["accession"] = df_filtered.saccver.str.strip()
    df_filtered["sample"] = sample_name
    df_filtered["sample"] = pd.Categorical(df_filtered["sample"])
    df_filtered["sample_segment"] = df_filtered.qaccver.str.extract(r".+_(\d)$").astype(
        "category"
    )
    df_filtered["sample_segment"] = pd.Categorical(df_filtered["sample_segment"])
    segments = df_filtered["sample_segment"].unique()
    df_filtered["subtype_from_match_title"] = (
        df_filtered["stitle"].str.extract(regex_subtype_pattern).astype("category")
    )
    df_filtered["subtype_from_match_title"] = df_filtered["subtype_from_match_title"]
    logging.info(
        f"{sample_name} | Merging NCBI Influenza DB genome metadata with BLAST results on accession."
    )
    df_merge = pd.merge(df_filtered, df_metadata, on="accession", how="left")
    del df_filtered
    del df_metadata
    df_merge["subtype"] = df_merge["subtype"].combine_first(df_merge["subtype_from_match_title"])
    df_merge = df_merge.sort_values(
        by=["sample_segment", "bitscore"], ascending=[True, False]
    ).set_index("sample_segment")
    subtype_results_summary = {}
    H_results = None
    N_results = None
    if "4" in df_merge.index:
        H_results = find_h_or_n_type(df_merge, "4")
        subtype_results_summary.update(H_results)
    if "6" in df_merge.index:
        N_results = find_h_or_n_type(df_merge, "6")
        subtype_results_summary.update(N_results)

    subtype_results_summary["sample"] = sample_name
    subtype_results_summary["subtype"] = get_subtype_value(H_results, N_results)
    dfs = []
    segments = df_merge.index.unique()
    for seg in segments:
        dfs.append(df_merge.loc[seg, :].head(top).reset_index())
    df_top_seg_matches = pd.concat(dfs)
    cols = pd.Series([x for x, _ in blast_results_report_columns])
    cols = cols[cols.isin(df_top_seg_matches.columns)]
    df_top_seg_matches = df_top_seg_matches[cols]
    return df_top_seg_matches, subtype_results_summary


def get_subtype_value(H_results: Optional[Dict], N_results: Optional[Dict]) -> str:
    subtype = ""
    if H_results is None and N_results is None:
        subtype = "-"
    elif H_results is not None and N_results is None:
        H: str = H_results.get("H_type", "")
        subtype = f"H{H}" if H != "" else "-"
    elif H_results is None and N_results is not None:
        N: str = N_results.get("N_type", "")
        subtype = f"N{N}" if N != "" else "-"
    else:
        H: str = H_results.get("H_type", "")
        N: str = N_results.get("N_type", "")
        if H == "" and N == "":
            subtype = "-"
        else:
            if H != "":
                H = f"H{H}"
            if N != "":
                N = f"N{N}"
            subtype = f"{H}{N}"
    return subtype


def find_h_or_n_type(df_merge, seg):
    assert seg in [
        "4",
        "6",
    ], "Can only determine H or N type from segments 4 or 6, respectively!"
    type_name = "H_type" if seg == "4" else "N_type"
    h_or_n = type_name[0]
    df_segment = df_merge.loc[seg, :]
    type_counts = df_segment.subtype.value_counts()
    df_type_counts = type_counts.index.str.extract(h_or_n + r"(\d+)")
    df_type_counts.columns = [type_name]
    df_type_counts["count"] = type_counts.values
    df_type_counts["subtype"] = type_counts.index
    df_type_counts = df_type_counts[~pd.isnull(df_type_counts[type_name])]
    logging.debug(f"{df_type_counts}")
    type_to_count = defaultdict(int)
    for _, x in df_type_counts.iterrows():
        if pd.isna(x[type_name]):
            continue
        type_to_count[x[type_name]] += x["count"]
    type_to_count = [(h, c) for h, c in type_to_count.items()]
    type_to_count.sort(key=lambda x: x[1], reverse=True)
    top_type, top_type_count = type_to_count[0]
    total_count = type_counts.sum()
    logging.info(
        f"{h_or_n}{top_type} n={top_type_count}/{total_count} ({top_type_count / total_count:.1%})"
    )
    type_mask = df_segment.subtype.str.match(
        r".*" + h_or_n + top_type + r".*", na=False
    )
    type_mask[pd.isnull(type_mask)] = False
    df_seg_top_type = df_segment[type_mask]
    top_result: pd.Series = [r for _, r in df_seg_top_type.head(1).iterrows()][0]
    results_summary = {
        f"{h_or_n}_type": top_type,
        f"{h_or_n}_sample_segment_length": top_result["qlen"],
        f"{h_or_n}_top_pident": top_result["pident"],
        f"{h_or_n}_top_mismatch": top_result["mismatch"],
        f"{h_or_n}_top_gaps": top_result["gapopen"],
        f"{h_or_n}_top_bitscore": top_result["bitscore"],
        f"{h_or_n}_top_align_length": top_result["length"],
        f"{h_or_n}_top_accession": top_result["accession"],
        f"{h_or_n}_top_host": top_result["host"],
        f"{h_or_n}_top_country": top_result["country"],
        f"{h_or_n}_top_date": top_result["date"],
        f"{h_or_n}_top_seq_length": top_result["slen"],
        f"{h_or_n}_virus_name": top_result["virus_name"],
        f"{h_or_n}_NCBI_Influenza_DB_subtype_matches": top_type_count,
        f"{h_or_n}_NCBI_Influenza_DB_total_matches": total_count,
        f"{h_or_n}_NCBI_Influenza_DB_proportion_matches": top_type_count / total_count,
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
        ("accession", str),
        ("host", "category"),
        ("segment", "category"),
        ("subtype", "str"),
        ("country", "category"),
        ("date", "category"),
        ("seq_length", "uint16"),
        ("virus_name", "category"),
        ("age", "category"),
        ("gender", "category"),
        ("group_id", "category"),
    ]
    df_md = pd.read_csv(
        flu_metadata,
        sep="\t",
        names=[name for name, _ in md_cols],
        dtype={name: t for name, t in md_cols},
    )
    unique_subtypes = df_md.subtype.unique()
    unique_subtypes = unique_subtypes[~pd.isna(unique_subtypes)]
    logging.info(
        f"Parsed Influenza metadata file into DataFrame with n={df_md.shape[0]} rows and n={df_md.shape[1]} columns. There are {unique_subtypes.size} unique subtypes. "
    )
    regex_subtype_pattern = r"\((H\d+N\d+|" + "|".join(list(unique_subtypes)) + r")\)"

    if threads > 1:
        pool = Pool(processes=threads)
        logging.info(
            f"Initialized multiprocessing pool with {threads} processes. Submitting async parsing jobs."
        )
        async_objects = [
            pool.apply_async(
                parse_blast_result,
                (blast_result, df_md, regex_subtype_pattern),
                dict(top=top, pident_threshold=pident_threshold),
            )
            for blast_result in blast_results
        ]
        logging.info(f"Getting async results...")
        results = [x.get() for x in async_objects]
        logging.info(
            f'Got {len(results)} async parsing results. Merging into report "{excel_report}".'
        )
    else:
        results = [
            parse_blast_result(blast_result, df_md, regex_subtype_pattern, top=top, pident_threshold=pident_threshold,
                               min_aln_length=min_aln_length) for blast_result in blast_results]
    dfs_blast = []
    all_subtype_results = {}
    for parsed_result in results:
        if parsed_result is None:
            continue
        df_blast, subtype_results_summary = parsed_result
        if df_blast is not None:
            dfs_blast.append(df_blast)
        sample = subtype_results_summary["sample"]
        all_subtype_results[sample] = subtype_results_summary
    df_subtype_results = pd.DataFrame(all_subtype_results).transpose()
    df_all_blast = pd.concat(dfs_blast).rename(
        columns={k: v for k, v in blast_results_report_columns}
    )
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
    if not get_top_ref:
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
        df_ref_id = df_all_blast[
            ['Sample', 'Sample Genome Segment Number', 'Reference NCBI Accession', 'BLASTN Bitscore',
             'Reference Sequence ID']]
        df_ref_id = df_ref_id.reset_index(drop=True)
        df_ref_id.loc[df_ref_id['Reference NCBI Accession'].isna(), 'Reference NCBI Accession'] = df_ref_id[
            'Reference Sequence ID']
        df_ref_id['Reference NCBI Accession'] = df_ref_id['Reference NCBI Accession'].str.strip()
        df_ref_id['Sample Genome Segment Number'] = df_ref_id['Sample Genome Segment Number']. \
            apply(lambda x: influenza_segment[int(x)])
        df_ref_id.to_csv(sample_name + ".topsegments.csv", header=True, index=False)


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
