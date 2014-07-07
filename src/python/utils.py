import ConfigParser
from bed_line import BedLine
import numpy as np
import csv

codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

def is_valid_nuc(nuc):
    valid_nucs = ['A', 'C', 'T', 'G', 'N']
    is_valid = nuc in valid_nucs
    return is_valid


def get_config(section):
    """Returns the config object"""
    cfg = ConfigParser.ConfigParser()
    cfg.read('config.cfg')
    cfg_options = dict(cfg.items(section))
    return cfg_options


def bed_generator(bed_path):
    with open(bed_path) as handle:
        bed_reader = csv.reader(handle, delimiter='\t')
        for line in bed_reader:
            yield BedLine(line)


def codon2aa(codon):
    codon = codon.upper()  # convert to upper case
    aa = codon_table[codon]
    return aa


def translate_seq(seq):
    seq_len = len(seq)
    if seq_len % 3:
        # coding sequence is not a multiple of three
        raise ValueError('Coding sequences should of multiple of three')

    protein_aa = [codon2aa(seq[i:i+3]) for i in range(0, seq_len, 3)]
    protein_seq = ''.join(protein_aa)
    return protein_seq


def lzip(*args):
    result = [[] for i in range(len(args[0]))]
    for arg in args:
        for i, ele in enumerate(arg):
            result[i].append(ele)
    return result


def cummin(x):
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(n)
    i = np.arange(1, n+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


