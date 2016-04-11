#!/usr/bin/env python
""" This script fetches sequence from an indexed FASTA for genes specified
in the provided BED file.
"""
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

import prob2020.python.utils as utils
import prob2020.python.gene_sequence as gs

# actually important imports
import pysam
import argparse
import logging
import datetime

logger = logging.getLogger(__name__)  # module logger


def parse_arguments():
    info = 'Extracts gene sequences from a genomic FASTA file'
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
                        default='',
                        help='Path to log file. (accepts stdout)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

    # program arguments
    help_str = 'Human genome FASTA file'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'BED file annotation of genes'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Output a single FASTA file with gene sequences'
    parser.add_argument('-o', '--output',
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

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))

    return vars(args)


def main(opts):
    # read bed file, extract gene sequence from genome, write to fasta
    genome_fa = pysam.Fastafile(opts['input'])
    with open(opts['output'], 'w') as handle:
        for bed_row in utils.bed_generator(opts['bed']):
            fasta_seq = gs.fetch_gene_fasta(bed_row, genome_fa)
            handle.write(fasta_seq)
    genome_fa.close()


def cli_main():
    opts = parse_arguments()
    main(opts)

if __name__ == "__main__":
    cli_main()
