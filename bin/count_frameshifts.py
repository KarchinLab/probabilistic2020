#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import permutation2020.python.utils as utils
import pandas as pd
import argparse


def keep_frameshifts(mut_df):
    # calculate length
    mut_df['indel len'] = mut_df['End_Position'] - mut_df['Start_Position']

    # filter out non-indels
    mut_df = mut_df[(mut_df['Reference_Allele']=='-') | (mut_df['Tumor_Allele']=='-')]

    # filter out non-frameshift indels
    mut_df = mut_df[(mut_df['indel len']%3)>0]

    return mut_df


def get_frameshift_lengths(num_bins):
    fs_len = []
    i = 1
    tmp_bins = 0
    while(tmp_bins<num_bins):
        if i%3:
            fs_len.append(i)
            tmp_bins += 1
        i += 1
    return fs_len


def count_frameshifts(mut_df,
                      bed_path,
                      num_bins,
                      use_unmapped=False):
    fs_cts = {}  # frameshift count information for each gene
    fs_df = keep_frameshifts(mut_df)
    fs_lens = get_frameshift_lengths(num_bins)

    for bed in utils.bed_generator(bed_path):
        gene_df = fs_df[fs_df['Gene']==bed.gene_name]

        # find it frameshift actually is on gene annotation
        fs_pos = []
        for ix, row in gene_df.iterrows():
            indel_pos = [row['Start_Position'], row['End_Position']]
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], indel_pos)
            fs_pos.append(coding_pos)

        # mark frameshifts that could not be mapped to reference tx
        gene_df['unmapped'] = [(1 if x is None else 0) for x in fs_pos]
        total_fs = len(gene_df)
        unmapped_fs = len(gene_df[gene_df['unmapped']==1])

        # filter out frameshifts that did not match reference tx
        if not use_unmapped:
            gene_df = gene_df[gene_df['unmapped']==0]

        # count frameshifts
        fs_len_cts = gene_df['indel len'].value_counts()[fs_lens]
        fs_len_cts = fs_len_cts.fillna(0).astype(int)

        fs_cts[bed.gene_name] = fs_len_cts.tolist() + [total_fs, unmapped_fs]

    fs_cts_df = pd.DataFrame.from_dict(fs_cts, orient='index')
    cols = map(str, fs_lens) + ['total', 'unmapped']
    fs_cts_df.columns = cols
    return fs_cts_df


def parse_arguments():
    info = 'Counts the number of frameshifts in each gene stratified by length'
    parser = argparse.ArgumentParser(description=info)
    help_str = 'Mutation file containing frameshift information'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Gene annotation in BED format with a single reference transcript'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Number of bins to categorize framshift lengths by'
    parser.add_argument('-bins', '--bins',
                        type=int, default=10,
                        help=help_str)
    help_str = 'Use frameshifts that could not be placed onto the reference transcript'
    parser.add_argument('-u', '--use-unmapped',
                        action='store_true', default=False,
                        help=help_str)
    help_str = 'Output of frameshift counts'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_csv(opts['mutations'], sep='\t')
    df['Start_Position'] = df['Start_Position'] - 1  # convert to 0-based coord
    fs = count_frameshifts(df, opts['bed'], opts['bins'], opts['use_unmapped'])
    fs.to_csv(opts['output'], sep='\t')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
