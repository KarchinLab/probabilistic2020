#!/usr/bin/env python
"""This scripts checks whether mutations are specified correctly.

Specifically, this script tests mutations to check whether they are
reported as being on the positive strand or on the coding strand.
"""
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# actually important imports
import prob2020.python.utils as utils
import pysam
import pandas as pd
import argparse
import logging
import datetime

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    info = 'Checks mutations to see what strand they are reported on and for unmapped mutations.'
    parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='stdout',
                        help='Path to log file. (Default: stdout)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

    # program arguments
    help_str = 'Human genome FASTA file'
    parser.add_argument('-f', '--fasta',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Text file specifying mutations in the format required for permutation test'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'BED file of reference transcripts'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Save mutations that could not be found on the reference transcript'
    parser.add_argument('-u', '--unmapped',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level,
                        verbose=args.verbose)  # start logging

    return vars(args)


def detect_coordinates(mut_df, genome_fa):
    # detect problems with using 0-based coordinates
    zero_len_count = 0
    num_snv = 0
    matching_ref = [0, 0]
    matching_pair = [0, 0]
    for ix, row in mut_df.iterrows():
        if (row['End_Position'] - row['Start_Position']) == 0:
            zero_len_count += 1
        no_shift_seq = genome_fa.fetch(reference=row['Chromosome'],
                                       start=row['Start_Position'],
                                       end=row['End_Position'])
        minus_1_seq = genome_fa.fetch(reference=row['Chromosome'],
                                      start=row['Start_Position']-1,
                                      end=row['End_Position'])
        seqs = [minus_1_seq, no_shift_seq]

        if len(row['Reference_Allele']) == 1 and row['Reference_Allele'] != '-':
            num_snv += 1

        for i in range(len(seqs)):
            if seqs[i].upper() == row['Reference_Allele'].upper() and len(row['Reference_Allele']) == 1:
                matching_ref[i] += 1
            elif seqs[i].upper() == utils.rev_comp(row['Reference_Allele']).upper() and len(row['Reference_Allele']) == 1:
                matching_pair[i] += 1

    # return coordinate type
    num_mut = len(mut_df)
    zero_len_pct = zero_len_count / float(num_mut)
    matching_pair_pct = map(lambda x: x / float(num_snv), matching_pair)
    matching_pct = map(lambda x: x / float(num_snv), matching_ref)
    logger.info('{0:.2f}%% for {1} tested mutations had zero length'.format(100*zero_len_pct, num_mut))
    logger.info('{0} for {1} did match the + strand reference genome'.format(matching_pct, num_snv))
    logger.info('{0} for {1} did match the - strand reference genome'.format(matching_pair_pct, num_snv))
    if zero_len_pct > .3:
        logger.info('1-based coordinate system likely used.')
        if matching_pair_pct[1] > .25:
            logger.info('Mutations likely reported on the genes\'s coding strand')
            return 1, 'coding'
        else:
            logger.info('Mutations likely reported on the genes\'s + strand')
            return 1, '+'
    elif (matching_ref[0] + matching_pair[0]) > (matching_ref[1] + matching_pair[1]):
        logger.info('0-based coordinate system likely used.')
        if matching_pair_pct[1] > .25:
            logger.info('Mutations likely reported on the genes\'s coding strand')
            return 0, 'coding'
        else:
            logger.info('Mutations likely reported on the genes\'s + strand')
            return 0, '+'
    else:
        logger.info('1-based coordinate system likely used.')
        if matching_pair_pct[1] > .25:
            logger.info('Mutations likely reported on the genes\'s coding strand')
            return 1, 'coding'
        else:
            logger.info('Mutations likely reported on the genes\'s + strand')
            return 1, '+'


def main(opts):
    # read in mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))
    mut_df = utils._fix_mutation_df(mut_df)

    # read genome fasta file
    genome_fa = pysam.Fastafile(opts['fasta'])

    # read BED file for transcripts
    bed_dict = utils.read_bed(opts['bed'], [])
    gene2bed = {item.gene_name: item
                for bed_list in bed_dict.values()
                for item in bed_list}

    # group mutations by gene
    mut_grpby = mut_df.groupby('Gene')
    unmapped_mut_list = []
    for i, mut_info in mut_grpby:
        gene_name = mut_info['Gene'].iloc[0]

        # try to find tx annotation for gene
        bed = None
        try:
            bed = gene2bed[gene_name]
        except KeyError:
            pass

        if bed:
            # get coding positions, mutations unmapped to the reference tx will have
            # NA for a coding position
            for ix, row in mut_info.iterrows():
                coding_pos = bed.query_position(bed.strand,
                                                row['Chromosome'],
                                                row['Start_Position'])
                if not coding_pos:
                    unmapped_mut_list.append(row.tolist())
        else:
            #unmapped_mut_df = pd.concat([unmapped_mut_df, mut_info])
            unmapped_mut_list += mut_info.values.tolist()

    # save the unmapped mutations to a file
    unmapped_mut_df = pd.DataFrame(unmapped_mut_list, columns=mut_df.columns)
    logger.info('{0} mutations were unmappable to a '
                'reference transcript'.format(len(unmapped_mut_df)))
    unmapped_mut_df.to_csv(opts['unmapped'], sep='\t', index=False)

    coord_base, coord_strand = detect_coordinates(mut_df, genome_fa)
    logger.info('RESULT: {0}-based coordinates, positions reported on {1} strand'.format(coord_base, coord_strand))

    genome_fa.close()


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
