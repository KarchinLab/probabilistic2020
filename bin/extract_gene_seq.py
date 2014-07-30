#!/usr/bin/env python
""" This script fetches sequence from an indexed FASTA for genes specified
in the provided BED file.
"""
# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../src/python'))

# actually important imports
import utils
import pysam
import argparse
import logging
import datetime

logger = logging.getLogger(__name__)  # module logger

def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """

    if not log_file:
        # create log directory if it doesn't exist
        file_dir = os.path.dirname(os.path.realpath(__file__))
        log_dir = os.path.join(file_dir, '../log/')
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def rev_comp(seq):
    """Get reverse complement of sequence.

    rev_comp will maintain the case of the sequence.

    Parameters
    ----------
    seq : str
        nucleotide sequence. valid {a, c, t, g, n}

    Returns
    -------
    rev_comp_seq : str
        reverse complement of sequence
    """
    base_pairing = {'A': 'T',
                    'T': 'A',
                    'a': 't',
                    't': 'a',
                    'C': 'G',
                    'G': 'C',
                    'c': 'g',
                    'g': 'c',
                    'n': 'n',
                    'N': 'N'}
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base_pairing[s] for s in rev_seq])
    return rev_comp_seq


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
    if strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end,
                             end=end+2)
        ss_seq = rev_comp(ss_seq)
    elif strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-2,
                             end=start)
    ss_fasta = '>{0};exon{1};5SS\n{2}\n'.format(gene_name,
                                                exon_num, ss_seq)
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

    if strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end,
                             end=end+2)
    elif strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-2,
                             end=start)
        ss_seq = rev_comp(ss_seq)
    ss_fasta = '>{0};exon{1};3SS\n{2}\n'.format(gene_name,
                                                exon_num, ss_seq)
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
            exon_seq = rev_comp(exon_seq)
        exon_fasta = '>{0};exon{1}\n{2}\n'.format(gene_bed.gene_name,
                                                  i, exon_seq)

        # get splice site sequence
        if len(exons) == 1:
            # splice sites don't matter if there is no splicing
            ss_fasta = ''
        elif i == 0:
            # first exon only, get 3' SS
            ss_fasta = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1]).upper()
        elif i == (len(exons) - 1):
            # last exon only, get 5' SS
            ss_fasta = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1]).upper()
        else:
            # middle exon, get bot 5' and 3' SS
            fasta_3ss = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1]).upper()
            fasta_5ss = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1]).upper()
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
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging

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
