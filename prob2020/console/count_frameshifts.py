#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import prob2020.python.utils as utils
import prob2020.python.indel as indel
import pandas as pd
import argparse


def count_frameshifts(mut_df,
                      bed_path,
                      num_bins,
                      num_samples,
                      use_unmapped=False):
    fs_cts = {}  # frameshift count information for each gene
    fs_df = indel.keep_frameshifts(mut_df)
    fs_lens = indel.get_frameshift_lengths(num_bins)

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
        all_len_counts = gene_df['indel len'].value_counts()
        fs_len_cts = all_len_counts[fs_lens]
        fs_len_cts[max(fs_len_cts.index)]= all_len_counts[all_len_counts.index>=fs_lens[-1]].sum()
        fs_len_cts = fs_len_cts.fillna(0).astype(int)

        # get length of gene
        gene_len = bed.cds_len
        gene_bases_at_risk = gene_len * num_samples

        info = [total_fs, unmapped_fs, gene_len, gene_bases_at_risk]
        fs_cts[bed.gene_name] = fs_len_cts.tolist() + info

    fs_cts_df = pd.DataFrame.from_dict(fs_cts, orient='index')
    cols = list(map(str, fs_lens)) + ['total', 'unmapped', 'gene length', 'bases at risk']
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
                        type=int, default=5,
                        help=help_str)
    help_str = 'Number of sequenced samples'
    parser.add_argument('-n', '--sample-number',
                        type=int, required=True,
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
    # read in data
    df = pd.read_csv(opts['mutations'], sep='\t')
    df['Start_Position'] = df['Start_Position'] - 1  # convert to 0-based coord

    # count frameshifts
    fs = count_frameshifts(df, opts['bed'], opts['bins'],
                           opts['sample_number'], opts['use_unmapped'])

    # save results
    if opts['output']:
        fs.to_csv(opts['output'], sep='\t')

    return fs


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
