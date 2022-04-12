#!/usr/bin/env python

import pandas as pd
import click
import matplotlib.pyplot as plt

NT_COLOURS = {
    'A': '#414cc1',
    'T': '#ead233',
    'G': '#3dc62d',
    'C': '#d3341f',
}


def parse_vcf(vcf_path) -> pd.DataFrame:
    df = pd.read_table(vcf_path,
                       comment='#',
                       header=None,
                       names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample'])
    return df


def read_depths(fpath) -> pd.DataFrame:
    df = pd.read_table(fpath,
                       names=['reference', 'pos', 'depth'],
                       header=None)
    return df


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
    return '; '.join([f'{xs[0]}-{xs[-1]}' for xs in coords])


def depth_plot(ax,
               df,
               low=3,
               sample_name='SAMPLE',
               segment='1',
               highlight_low_cov=True,
               highlight_no_cov=True) -> str:
    plt.sca(ax)
    genome = sample_name
    reference = df.reference.values[0]
    dflow = df.copy()
    low_depth = (dflow.depth < low)
    dflow.loc[low_depth, 'depth'] = df.depth.max()
    df0 = df.copy()
    zero_depth = (df0.depth == 0)
    df0.loc[zero_depth, 'depth'] = df.depth.max()
    ax.set_title(f'Sample "{genome}", segment {segment} mapped to reference "{reference}"\n')
    ax.set_ylabel('Depth')
    ax.set_xlabel('Position')
    ax.set_ylim(top=df.depth.max())
    ax.set_xlim(left=1, right=df.pos.max())
    ax.fill_between('pos', 'depth', 0, data=df, color='darkgrey')
    if highlight_low_cov:
        ax.fill_between('pos', 'depth', df.depth, where=dflow.depth > df.depth, color='yellow', data=dflow)
    if highlight_no_cov:
        ax.fill_between('pos', 'depth', df.depth, where=df0.depth > df.depth, color='red', data=df0)

    return (f'Mean (median) coverage: {df.depth.mean():.1f}X ({df.depth.median():.1f}X)\n'
            f'Genome coverage (>= {low}X): {(df.depth >= low).sum() / df.shape[0]:.1%}\n'
            f'Low coverage positions (<{low}X): {low_depth.sum():.0f}\n'
            f'No coverage positions (0X): {zero_depth.sum():.0f}\n'
            f'Low coverage regions (<{low}X): {get_interval_coords(df, low - 1)}\n'
            f'No coverage regions (0X): {get_interval_coords(df, 0)}\n')


@click.command()
@click.option('-d', '--depths-file', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-pdf', required=True, type=click.Path(exists=False))
@click.option('-v', '--vcf-file', type=click.Path(exists=True))
@click.option('-l', '--low-coverage', default=3, type=int)
@click.option('--log-scale-y', is_flag=True)
@click.option('-w', '--width', default=10, type=int)
@click.option('-h', '--height', default=5, type=int)
@click.option('--sample-name', default='SAMPLE', type=str)
@click.option('--segment', default='1', type=str)
@click.option('--no-highlight', is_flag=True)
def main(depths_file,
         output_pdf,
         vcf_file,
         low_coverage,
         log_scale_y,
         width,
         height,
         sample_name,
         segment,
         no_highlight):
    highlight_low_no_cov_regions = not no_highlight
    df = read_depths(depths_file)
    df_vcf = None
    if vcf_file:
        df_vcf = parse_vcf(vcf_file)
        print(df_vcf)
    fig, ax = plt.subplots(1, figsize=(width, height))
    if log_scale_y:
        from matplotlib.ticker import ScalarFormatter
        ax.set_yscale('log', nonpositive='clip')
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter)
    bottom_desc = depth_plot(ax,
                             df,
                             low=low_coverage,
                             sample_name=sample_name,
                             segment=segment,
                             highlight_low_cov=highlight_low_no_cov_regions,
                             highlight_no_cov=highlight_low_no_cov_regions)
    # add coverage stats description below plot
    fig.text(0.1, -0.15, bottom_desc, fontsize='small')
    print(bottom_desc)
    if df_vcf is not None:
        df = df.set_index('pos')
        for idx, row in df_vcf.iterrows():
            depth = df.loc[row.POS, 'depth']
            color = NT_COLOURS.get(row.ALT, '#555555')
            ax.axvline(x=row.POS, color=color, alpha=0.4)
            ax.text(row.POS,
                    depth,
                    f'{row.REF}{row.POS}{row.ALT}\nDP={depth}',
                    fontsize='xx-small')
    fig.savefig(output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
