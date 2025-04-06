#!/usr/bin/env python

import logging
from pathlib import Path

import polars as pl
import typer
from rich.logging import RichHandler

VERSION = "2025.04.0"

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


def version_callback(value: bool):
    if value:
        typer.echo(f"{VERSION}")
        raise typer.Exit()


@app.command()
def main(
        nextclade_tsv_input_dir: Path = typer.Argument(
            ...,
            exists=True,
            dir_okay=True,
            help="Directory with Nextclade TSV files."
        ),
        manifest_csv: Path = typer.Argument(
            ...,
            exists=True,
            help="CSV with sample name, dataset name and tag (if applicable) and Nextclade output TSV filename"
        ),
        nextclade_tsv_output: Path = typer.Argument(
            "nextclade.tsv",
            help="Output TSV of aggregated and filtered Nextclade results"
        ),
        verbose: bool = typer.Option(False, "--verbose", "-v", is_flag=True, help="Enable verbose logging"),
        version: bool = typer.Option(None, "--version", callback=version_callback, is_eager=True),
):
    """Aggregate positive results from multiple Nextclade analyses across datasets and samples into a single report.

    """
    init_logging(verbose)
    nextclade_tsv_filename2path = {p.name: p for p in nextclade_tsv_input_dir.glob("*.nextclade.tsv")}
    df_manifest = pl.read_csv(
        manifest_csv,
        has_header=False,
        new_columns=["sample", "dataset_name", "dataset_tag", "nextclade_tsv_filename"]
    )
    dfs = []
    for row in df_manifest.iter_rows(named=True):
        dataset_name = row["dataset_name"]
        dataset_tag = row["dataset_tag"]
        if dataset_tag == "null":
            dataset_tag = "latest"
        sample = row["sample"]
        nextclade_tsv_filename = row["nextclade_tsv_filename"]
        nextclade_tsv_path = nextclade_tsv_filename2path[nextclade_tsv_filename]
        df = pl.read_csv(nextclade_tsv_path, separator="\t", infer_schema_length=0)
        df = df.with_columns([
            pl.col(col).cast(pl.Utf8).alias(col) for col in df.columns
        ])
        df = df.with_columns([
            pl.lit(sample).alias("sample"),
            pl.lit(dataset_name).alias("dataset_name"),
            pl.lit(dataset_tag).alias("dataset_tag")
        ])
        logger.info(f"{sample}: {df.shape=}")
        dfs.append(df)
    all_cols = {col: df[col].dtype for df in dfs for col in df.columns}
    first_ordered_cols = [
        'sample',
        'seqName',
        'dataset_name',
        'dataset_tag',
        'clade',
        'subclade',
    ]
    rest_cols = [col for col in all_cols if col not in first_ordered_cols]
    rest_cols.sort()
    ordered_cols = first_ordered_cols + rest_cols
    logger.info(
        "Ensuring that all Nextclade TSV outputs have the same columns in the same order. Adding extra columns if necessary and reordering.")
    dfs_all_cols = []
    for df in dfs:
        for col, col_dtype in all_cols.items():
            if col not in df.columns:
                null_val = None
                if col_dtype == pl.Int64:
                    null_val = 0
                elif col_dtype == pl.Float64:
                    null_val = 0.0
                elif col_dtype == pl.Utf8:
                    null_val = ''
                df = df.with_columns(pl.lit(null_val).alias(col))
        df = df.select(ordered_cols)
        dfs_all_cols.append(df)
    logger.info(f"Concatenating {len(dfs_all_cols)} Nextclade reports into single report.")
    df_concat = pl.concat(dfs_all_cols)
    logger.info(f"Unfiltered Nextclade report size {df_concat.shape}")
    df_filtered = df_concat.filter(pl.col('errors').is_null() | ~(pl.col('errors').str.contains('Unable')))
    logger.info(f"Filtered Nextclade report size {df_filtered.shape}")
    df_filtered = df_filtered.sort(by=['sample', 'seqName', 'dataset_name'], descending=False)
    df_filtered.write_csv(nextclade_tsv_output, separator='\t')
    logger.info(f"Wrote aggregated and filtered Nextclade output to '{nextclade_tsv_output}'")


if __name__ == "__main__":
    app()
