#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

from prob2020.python.bed_line import BedLine
import pandas as pd
import argparse
import IPython


def parse_arguments():
    info = 'Only keeps the longest gene from the BED file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help='input bed file')
    parser.add_argument('-d', '--duplicated',
                        type=str, required=True,
                        help='output file for list of duplicated genes')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='output bed file keeping the longest tx')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_csv(opts['bed'], sep='\t', header=None)
    gene_lengths = []
    for i, row in df.iterrows():
        gene_lengths.append(BedLine(row.tolist()).cds_len)
    df['gene_length'] = gene_lengths
    gene_occurs = df[3].value_counts()
    duplicated_cts = gene_occurs[gene_occurs>1]
    dup_genes = duplicated_cts.index
    dup_df = df[df[3].isin(dup_genes)].copy()
    grp = dup_df.groupby(3)
    small_tx = grp.apply(lambda d: d.index[d.gene_length.argmin()])
    dedup_df = df[~df.index.isin(small_tx)]

    dedup_df.drop(['gene_length'], axis=1).to_csv(opts['output'], sep='\t', index=False, header=False)
    duplicated_cts.to_csv(opts['duplicated'], sep='\t')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
