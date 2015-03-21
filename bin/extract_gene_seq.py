#!/usr/bin/env python
""" This script fetches sequence from an indexed FASTA for genes specified
in the provided BED file.
"""
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import prob2020.python.utils as utils

# actually important imports
import pysam
import argparse
import logging
import datetime

logger = logging.getLogger(__name__)  # module logger

def _fetch_5ss_fasta(fasta, gene_name, exon_num,
                     chrom, strand, start, end):
    """Retreives the 5' SS sequence flanking the specified exon.

    Returns a string in fasta format with the first line containing
    a ">" and the second line contains the two base pairs of 5' SS.

    Parameters
    ----------
    fasta : pysam.Fastafile
        fasta object from pysam
    gene_name : str
        gene name used for fasta seq id
    exon_num : int
        the `exon_num` exon, used for seq id
    chrom : str
        chromsome
    strand : str
        strand, {'+', '-'}
    start : int
        0-based start position
    end : int
        0-based end position

    Returns
    -------
    ss_fasta : str
        string in fasta format with first line being seq id
    """
    if strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end-1,
                             end=end+3)
    elif strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-3,
                             end=start+1)
        ss_seq = utils.rev_comp(ss_seq)
    ss_fasta = '>{0};exon{1};5SS\n{2}\n'.format(gene_name,
                                                exon_num,
                                                ss_seq.upper())
    return ss_fasta


def _fetch_3ss_fasta(fasta, gene_name, exon_num,
                     chrom, strand, start, end):
    """Retreives the 3' SS sequence flanking the specified exon.

    Returns a string in fasta format with the first line containing
    a ">" and the second line contains the two base pairs of 3' SS.

    Parameters
    ----------
    fasta : pysam.Fastafile
        fasta object from pysam
    gene_name : str
        gene name used for fasta seq id
    exon_num : int
        the `exon_num` exon, used for seq id
    chrom : str
        chromsome
    strand : str
        strand, {'+', '-'}
    start : int
        0-based start position
    end : int
        0-based end position

    Returns
    -------
    ss_fasta : str
        string in fasta format with first line being seq id
    """

    if strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end-1,
                             end=end+3)
        ss_seq = utils.rev_comp(ss_seq)
    elif strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-3,
                             end=start+1)
    ss_fasta = '>{0};exon{1};3SS\n{2}\n'.format(gene_name,
                                                exon_num,
                                                ss_seq.upper())
    return ss_fasta


def fetch_gene_fasta(gene_bed, fasta_obj):
    """Retreive gene sequences in FASTA format.

    Parameters
    ----------
    gene_bed : BedLine
        BedLine object representing a single gene
    fasta_obj : pysam.Fastafile
        fasta object for index retreival of sequence

    Returns
    -------
    gene_fasta : str
        sequence of gene in FASTA format
    """
    gene_fasta = ''
    strand = gene_bed.strand
    exons = gene_bed.get_exons()
    if strand == '-':
        exons.reverse()  # order exons 5' to 3', so reverse if '-' strand

    # iterate over exons
    for i, exon in enumerate(exons):
        exon_seq = fasta_obj.fetch(reference=gene_bed.chrom,
                                   start=exon[0],
                                   end=exon[1]).upper()
        if strand == '-':
            exon_seq = utils.rev_comp(exon_seq)
        exon_fasta = '>{0};exon{1}\n{2}\n'.format(gene_bed.gene_name,
                                                  i, exon_seq)

        # get splice site sequence
        if len(exons) == 1:
            # splice sites don't matter if there is no splicing
            ss_fasta = ''
        elif i == 0:
            # first exon only, get 3' SS
            ss_fasta = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1])
        elif i == (len(exons) - 1):
            # last exon only, get 5' SS
            ss_fasta = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1])
        else:
            # middle exon, get bot 5' and 3' SS
            fasta_3ss = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1])
            fasta_5ss = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1])
            ss_fasta = fasta_5ss + fasta_3ss

        gene_fasta += exon_fasta + ss_fasta

    return gene_fasta


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
            fasta_seq = fetch_gene_fasta(bed_row, genome_fa)
            handle.write(fasta_seq)
    genome_fa.close()


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
