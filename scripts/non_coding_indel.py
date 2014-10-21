import pandas as pd
import pysam
import csv
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        required=True, type=str,
                        help='File of indel mutations')
    parser.add_argument('-b', '--blacklist',
                        required=True, type=str,
                        help='List of genomic coordinates in BED format '
                             'that are not allowed for non-coding background')
    parser.add_argument('-bo', '--blacklist-output',
                        required=True, type=str,
                        help='Indel mutations that are in black list')
    parser.add_argument('-o', '--output',
                        required=True, type=str,
                        help='Indel mutations not in black list')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in INDEL mutations
    indels = pd.read_csv(opts['input'], sep='\t')

    # pysam tabix uses 1-based coordinates
    pysam.tabix_index(opts['blacklist'], force=True,
                      seq_col=0, start_col=1, end_col=2)

    # query black list to find INDELs with no hits
    non_coding_ixs, coding_ixs = [], []
    black_list = pysam.Tabixfile(opts['blacklist'])
    for i, row in indels.iterrows():
        result = black_list.fetch(reference=row['Chromosome'],
                                  start=row['Start_Position'],
                                  end=row['End_Position'])
        if not list(result):
            non_coding_ixs.append(i)
        else:
            coding_ixs.append(i)
    black_list.close()

    # save non-coding indels
    indels.ix[non_coding_ixs, :].to_csv(opts['output'], sep='\t', index=False)
    indels.ix[coding_ixs, :].to_csv(opts['blacklist_output'], sep='\t', index=False)


if __name__=="__main__":
    opts = parse_arguments()
    main(opts)
