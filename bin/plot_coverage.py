#!/usr/bin/env python

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import typer
from rich.logging import RichHandler

NT_COLOURS = {
    "A": "#414cc1",
    "T": "#ead233",
    "G": "#3dc62d",
    "C": "#d3341f",
}

SEGMENT_NAMES = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "M",
    "8": "NS",
}


def parse_vcf(vcf_path: Path) -> pd.DataFrame:
    df = pd.read_table(
        vcf_path,
        comment="#",
        header=None,
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "sample",
        ],
    )
    return df


def parse_depths(path: Path) -> pd.DataFrame:
    df_bed = read_mosdepth_perbase_bedgz(path)
    return expand_bed(df_bed)


def expand_bed(df_bed: pd.DataFrame) -> pd.DataFrame:
    arr = depth_array(df_bed)
    return pd.DataFrame(dict(pos=range(1, arr.size + 1), depth=arr))


def read_mosdepth_perbase_bedgz(path: Path) -> pd.DataFrame:
    df = pd.read_table(
        path, comment="#", header=None, names=["genome", "start_idx", "end_idx", "depth"]
    )
    return df


def depth_array(df: pd.DataFrame) -> np.ndarray:
    """Mosdepth per-base BED to depths array

    Data type of array is unsigned 32-bit int (`np.uint32`) (max value = 4294967295)

    >>> depth_array(pd.DataFrame(dict(start_idx=[0, 5, 10], end_idx=[5, 10, 15], depth=[1, 3, 5])))
    array([1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5], dtype=uint32)
    """
    arr = np.zeros(df.end_idx.max(), dtype=np.uint32)
    for row in df.itertuples():
        arr[row.start_idx: row.end_idx] = int(row.depth)
    return arr


def get_interval_coords(df, threshold=0):
    pos = df[df.depth <= threshold].pos
    coords = []
    for i, x in enumerate(pos):
        if coords:
            last = coords[-1][-1]
            if x == last + 1:
                coords[-1].append(x)
            else:
                coords.append([x])
        else:
            coords.append([x])
    return "; ".join([f"{xs[0]}-{xs[-1]}" for xs in coords])


def depth_plot(
        ax,
        df: pd.DataFrame,
        ref_id: str,
        low=3,
        sample_name="SAMPLE",
        segment="1",
        highlight_low_cov=True,
        highlight_no_cov=True,
) -> str:
    plt.sca(ax)
    dflow = df.copy()
    low_depth = dflow.depth < low
    dflow.loc[low_depth, "depth"] = df.depth.max()
    df0 = df.copy()
    zero_depth = df0.depth == 0
    df0.loc[zero_depth, "depth"] = df.depth.max()
    ax.set_title(
        f"Sample {sample_name} - Segment {segment} ({SEGMENT_NAMES[segment]})\nReference Sequence ID: {ref_id}"
    )
    ax.set_ylabel("Depth")
    ax.set_xlabel(f"Position in {ref_id}")
    ax.set_ylim(top=df.depth.max())
    ax.set_xlim(left=1, right=df.pos.max())
    ax.fill_between("pos", "depth", 0, data=df, color="darkgrey")
    if highlight_low_cov:
        ax.fill_between(
            "pos",
            "depth",
            df.depth,
            where=dflow.depth > df.depth,
            color="yellow",
            data=dflow,
        )
    if highlight_no_cov:
        ax.fill_between(
            "pos", "depth", df.depth, where=df0.depth > df.depth, color="red", data=df0
        )

    return (
        f"Mean (median) coverage: {df.depth.mean():.1f}X ({df.depth.median():.1f}X)\n"
        f"Genome coverage (>= {low}X): {(df.depth >= low).sum() / df.shape[0]:.1%}\n"
        f"Low coverage positions (<{low}X): {low_depth.sum():.0f}\n"
        f"No coverage positions (0X): {zero_depth.sum():.0f}\n"
        f"Low coverage regions (<{low}X): {get_interval_coords(df, low - 1)}\n"
        f"No coverage regions (0X): {get_interval_coords(df, 0)}\n"
    )


def init_logging():
    from rich.traceback import install

    install(show_locals=True, width=200, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def main(
        depths_file: Path = typer.Option(..., "-d", "--depths"),
        output_pdf: Path = typer.Option(..., "-o", "--output-pdf"),
        vcf_file: Path = typer.Option(..., "-v", "--vcf-file"),
        low_coverage: int = typer.Option(9),
        log_scale_y: bool = typer.Option(False),
        width: float = typer.Option(10.0),
        height: float = typer.Option(5.0),
        sample_name: str = typer.Option("SAMPLE"),
        segment: str = typer.Option("1"),
        ref_id: str = typer.Option(None),
        no_highlight: bool = typer.Option(False),
):
    init_logging()
    if "_" in segment:
        segment, _ = segment.split('_', maxsplit=1)
    highlight_low_no_cov_regions = not no_highlight
    df = parse_depths(depths_file)
    df_vcf = None
    if vcf_file:
        df_vcf = parse_vcf(vcf_file)
    fig, ax = plt.subplots(1, figsize=(width, height))
    if len(df):
        if log_scale_y:
            from matplotlib.ticker import ScalarFormatter

            ax.set_yscale("log", nonpositive="clip")
            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            ax.yaxis.set_major_formatter(formatter)
        bottom_desc = depth_plot(
            ax,
            df,
            ref_id=ref_id,
            low=low_coverage,
            sample_name=sample_name,
            segment=segment,
            highlight_low_cov=highlight_low_no_cov_regions,
            highlight_no_cov=highlight_low_no_cov_regions,
        )
        # add coverage stats description below plot
        fig.text(0.1, -0.15, bottom_desc, fontsize="small")
    if df_vcf is not None:
        df = df.set_index("pos")
        for idx, row in df_vcf.iterrows():
            depth = df.loc[row.POS, "depth"]
            color = NT_COLOURS.get(row.ALT, "#555555")
            ax.axvline(x=row.POS, color=color, alpha=0.4)
            ax.text(
                row.POS,
                depth,
                f"{row.REF}{row.POS}{row.ALT}\nDP={depth}",
                fontsize="xx-small",
            )
    if len(df) == 0:
        ax.set_title(f"No mapping (Empty Plot)\n")
    fig.savefig(output_pdf, bbox_inches="tight")


if __name__ == "__main__":
    typer.run(main)
