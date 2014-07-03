#!/usr/bin/env python
"""Create a BED file (using gene names) that describes the
gene's longest transcript as defined by SNVBox.

Transcripts with multiple locations are removed.
"""
import pandas as pd
import argparse


def parse_arguments():
    d = 'Returns a file with only one tx per gene and fixes tx name issues.'
    parser = argparse.ArgumentParser(description=d)
    help_str = ('BED output from UCSC table browser for transcripts '
                'identified by longest_snvbox_tx.sql script')
    parser.add_argument('-b', '--bed',
                        required=True, type=str,
                        help=help_str)
    help_str = 'Gene annotation output from unique_tx.py script'
    parser.add_argument('-g', '--genes',
                        required=True, type=str,
                        help=help_str)
    help_str = 'Path to save BED file describing genes rather than transcripts'
    parser.add_argument('-o', '--output',
                        required=True, type=str,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in table browser output
    bed_cols = ['chrom', 'chromStart', 'chromEnd', 'name',
                'score', 'strand', 'thickStart', 'thickEnd',
                'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    tx_df = pd.read_csv(opts['bed'], sep='\t',
                        header=None, names=bed_cols)

    # drop all rows that are part of duplicate positions (eg tx on chrX and
    # chrY)
    tx_grp = tx_df.groupby('name')
    tx_df = tx_df.set_index('name')  # set tx ids as indices
    tx_df = tx_df[tx_grp.apply(lambda df: len(df)==1)]  # drop rows if tx in multiple positions

    # read in genes defined by SNVBox
    gene_df = pd.read_csv(opts['genes'], sep='\t')

    # match gene names with tx in bed file by doing a table merge (like SQL)
    merged_df = pd.merge(gene_df, tx_df,
                         how='inner',
                         left_on='Transcript', right_index=True)
    merged_df = merged_df.rename(columns={'Gene Symbol': 'name'})

    # write BED file with the 'name' column now containing gene names rather
    # than transcript names
    merged_df[bed_cols].to_csv(opts['output'], sep='\t',
                               header=None, index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
