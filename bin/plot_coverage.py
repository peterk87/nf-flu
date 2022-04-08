#!/usr/bin/env python

import pandas as pd
import click
import matplotlib.style as mplstyle
import matplotlib.pyplot as plt


NT_COLOURS = {
	'A': '#414cc1',
	'T': '#ead233',
	'G': '#3dc62d',
	'C': '#d3341f',
}


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


@click.command()
@click.option('-d', '--depths-file', required=True, type=click.Path(exists=True))
@click.option('-o', '--output-pdf', required=True, type=click.Path(exists=False))
@click.option('-v', '--vcf-file', type=click.Path(exists=True))
@click.option('-l', '--low-coverage', default=3, type=int)
@click.option('--log-scale-y', is_flag=True)
@click.option('-w', '--width', default=12, type=int)
@click.option('-h', '--height', default=10, type=int)
def main(depths_file, 
		 output_pdf,
		 vcf_file,
		 low_coverage,
		 log_scale_y,
		 width,
		 height):
	
	mplstyle.use(['seaborn',])
	df = read_depths(depths_file)
	df_vcf = None
	if vcf_file:
		df_vcf = parse_vcf(vcf_file)
		print(df_vcf)
	fig, ax = plt.subplots(1, figsize=(width, height))
	if log_scale_y:
		from matplotlib.ticker import ScalarFormatter
		ax.set_yscale('log', nonposy='clip')
		formatter = ScalarFormatter()
		formatter.set_scientific(False)
		ax.yaxis.set_major_formatter(formatter)
	depth_plot(ax, df, low=low_coverage)
	if df_vcf is not None:
		df = df.set_index('pos')
		for idx, row in df_vcf.iterrows():
			depth = df.loc[row.POS, 'depth']
			color = NT_COLOURS.get(row.ALT, '#555555')
			ax.axvline(x=row.POS, color=color, alpha=0.5)
			ax.text(row.POS, depth, f'{row.REF}->{row.ALT}\nPostion: {row.POS}', fontsize='xx-small')
	fig.savefig(output_pdf, bbox_inches='tight')


if __name__ == '__main__':
	main()