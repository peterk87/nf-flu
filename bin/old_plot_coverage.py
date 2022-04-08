#!/usr/bin/env python

import pandas as pd
import matplotlib.style as mplstyle
import matplotlib.pyplot as plt
import argparse

NT_COLOURS = {
    'A': '#414cc1',
    'T': '#ead233',
    'G': '#3dc62d',
    'C': '#d3341f',
}


def parse_args(args=None):
    Description = "Plot Coverage Plot"
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument(
        "-d",
        "--depths-file",
        help="Coverage Depth Files",
    )
    parser.add_argument(
        "-o",
        "--output-pdf",
        help="Output Files.",
    )
    parser.add_argument(
        "-v",
        "--vcf-file",
        help="VCF files.",
    )
    parser.add_argument(
        "-l",
        "--low-coverage",
        default=3, type=int,
    )
    parser.add_argument(
        "-s",
        "--log-scale-y",
        default=True, type=bool,
    )
    return parser.parse_args(args)


def parse_vcf(vcf_path) -> pd.DataFrame:
    df = pd.read_table(vcf_path, comment='#', header=None,
                       names='CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample'.split())
    return df


def read_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names='genome reference pos depth'.split(),
                       header=None)
    return df


def depth_plot(ax, df, low=3):
    plt.sca(ax)
    genome = df.genome.values[0]
    reference = df.reference.values[0]
    dflow = df.copy()
    low_depth = (dflow.depth <= low)
    dflow.loc[low_depth, 'depth'] = df.depth.max()
    df0 = df.copy()
    zero_depth = (df0.depth == 0)
    df0.loc[zero_depth, 'depth'] = df.depth.max()
    ax.set_title(
        f'Sample "{genome}" mapped to reference "{reference}"\n'
        f'Mean coverage: {df.depth.mean():.1f}X\n'
        f'Low coverage positions (<=3X): {low_depth.sum():.0f}\n'
        f'Zero coverage positions (0X): {zero_depth.sum():.0f}'
    )
    ax.set_ylabel('Depth')
    ax.set_xlabel('Position')
    ax.set_ylim(top=df.depth.max())
    ax.set_xlim(left=1, right=df.pos.max())
    ax.fill_between('pos', 'depth', 0, data=df, color='darkgrey')
    ax.fill_between('pos', 'depth', df.depth, where=dflow.depth > df.depth, color='yellow', data=dflow)
    ax.fill_between('pos', 'depth', df.depth, where=df0.depth > df.depth, color='red', data=df0)


def main(args=None):
    args = parse_args(args)
    mplstyle.use(['seaborn', ])
    df = read_depths(args.depths_file)
    df_vcf = None
    if args.vcf_file:
        df_vcf = parse_vcf(args.vcf_file)
        print(df_vcf)
    fig, ax = plt.subplots(1, figsize=(12,10))
    if args.log_scale_y:
        from matplotlib.ticker import ScalarFormatter
        ax.set_yscale('log', nonposy='clip')
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter)
    depth_plot(ax, df, low=args.low_coverage)
    if df_vcf is not None:
        df = df.set_index('pos')
        for idx, row in df_vcf.iterrows():
            depth = df.loc[row.POS, 'depth']
            color = NT_COLOURS.get(row.ALT, '#555555')
            ax.axvline(x=row.POS, color=color, alpha=0.5)
            ax.text(row.POS, depth, f'{row.REF}->{row.ALT}\nPostion: {row.POS}', fontsize='xx-small')
    fig.savefig(args.output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
