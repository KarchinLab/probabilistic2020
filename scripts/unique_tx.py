#!/usr/bin/env python
"""The longest_snvbox_tx SQL script returns all transcripts
that are of maximum length. This script grabs the first occurrence
for a gene so that each gene has a single tx.
"""
import pandas as pd
import re
import argparse


def fix_tx_name(tx):
    """Remove the .[0-9] at the end of RefSeq IDs."""
    if tx.startswith('NM'):
        return re.sub('\.[0-9]$', '', tx)
    else:
        return tx


def parse_arguments():
    d = 'Returns a file with only one tx per gene and fixes tx name issues.'
    parser = argparse.ArgumentParser(description=d)
    help_str = 'Text output from longest_snvbox_tx.sql'
    parser.add_argument('-i', '--input',
                        required=True, type=str,
                        help=help_str)
    help_str = 'Path to save output'
    parser.add_argument('-o', '--output',
                        required=True, type=str,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # Fix refseq ids, then select first tx for every gene, finally save output
    df = pd.read_csv(opts['input'], sep='\t')
    df['Transcript'] = df['Transcript'].apply(fix_tx_name)
    one_tx = df.groupby('Gene Symbol').first()  # select first tx for a gene
    one_tx.to_csv(opts['output'], sep='\t')


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
