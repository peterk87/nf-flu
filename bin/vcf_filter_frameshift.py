#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from io import TextIOWrapper
from os import PathLike
from pathlib import Path
from typing import TextIO, Union, Tuple

import pandas as pd
import typer


VCF_COL_DTYPES: dict = dict(CHROM='category',
                            POS='uint32',
                            ID='category',
                            REF='category',
                            ALT='category',
                            QUAL=float,
                            FILTER='category',
                            INFO=str,
                            FORMAT=str)


def read_vcf(vcf_file: Union[PathLike, TextIOWrapper, TextIO]) -> Tuple[str, pd.DataFrame]:
    """Read VCF file into a DataFrame"""
    header = ''
    with vcf_file if isinstance(vcf_file, (TextIOWrapper, TextIO)) else open(vcf_file) as fh:
        vcf_cols = []
        for line in fh:
            header += line
            if line.startswith('#CHROM'):
                vcf_cols = line[1:].strip().split('\t')
                break
        df = pd.read_table(fh,
                           comment='#',
                           header=None,
                           names=vcf_cols,
                           dtype=VCF_COL_DTYPES)
    return header, df


def write_vcf(vcf_file: Union[PathLike, TextIOWrapper, TextIO], header: str, df: pd.DataFrame) -> None:
    with vcf_file if isinstance(vcf_file, (TextIOWrapper, TextIO)) else open(vcf_file, 'w') as fh:
        fh.write(header)
        df.to_csv(fh, sep='\t', header=False, index=False)


def main(input_vcf: Path = typer.Argument(None, help='VCF file to filter'),
         output_vcf: Path = typer.Argument(None, help='Ouput VCF file to filter'),
         verbose: bool = typer.Option(True)):
    """Filter variants that lead to frameshift mutations from a VCF file

    Any indels resulting in a frameshift (indel length not divisible by 3 (AA codon length)) are
    removed from the output VCF file.
    """
    if input_vcf and not input_vcf.is_file():
        log.warning(f'input_vcf not a file or stdin stream!')
        sys.exit(1)

    # reader = vcfpy.Reader.from_stream(input_vcf) if isinstance(input_vcf, TextIOWrapper) else vcfpy.Reader.from_path(input_vcf)
    header, df = read_vcf(input_vcf if input_vcf else sys.stdin)
    # writer = vcfpy.Writer.from_stream(output_vcf, reader.header) if isinstance(output_vcf, TextIOWrapper) else vcfpy.Writer.from_path(output_vcf, reader.header)
    # for idx, rec in reader.iterrows():
    potential_frameshift_mask = df.apply(lambda x:  (len(x.ALT) - len(x.REF)) % 3 == 0, axis=1)
    df_filtered = df[potential_frameshift_mask]

    # log any potential frameshift variants that are going to be filtered out in output VCF
    if (~potential_frameshift_mask).any():
        #log.info(f'{(~potential_frameshift_mask).sum()} frameshift variants filtered out')
        for _, row in df[~potential_frameshift_mask].iterrows():
            print(row)

    write_vcf(output_vcf if output_vcf else sys.stdout, header, df_filtered)


if __name__ == '__main__':
    typer.run(main)
