# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../cython'))

# normal imports
import ConfigParser
from bed_line import BedLine
import numpy as np
import csv
import itertools as it
import cutils


# global dictionary mapping codons to AA
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
    """Iterates through a BED file yielding parsed BED lines.

    Parameters
    ----------
    bed_path : str
        path to BED file

    Yields
    ------
    BedLine(line) : BedLine
        A BedLine object which has parsed the individual line in
        a BED file.
    """
    with open(bed_path) as handle:
        bed_reader = csv.reader(handle, delimiter='\t')
        for line in bed_reader:
            yield BedLine(line)


def codon2aa(codon):
    """Gets corresponding AA for a codon.

    Handles lower case as well as upper case.
    """
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
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.

    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.

    Parameters
    ----------
    pval : list or array
        list/array of p-values

    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
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


def get_aa_mut_info(coding_pos, somatic_base, gene_seq):
    """Retrieves relevant information about the effect of a somatic
    SNV on the amino acid of a gene.

    Information includes the germline codon, somatic codon, codon
    position, germline AA, and somatic AA.

    Parameters
    ----------
    coding_pos : iterable of ints
        Contains the base position (0-based) of the mutations
    somatic_base : list of str
        Contains the somatic nucleotide for the mutations
    gene_seq : GeneSequence
        gene sequence

    Returns
    -------
    aa_info : dict
        information about the somatic mutation effect on AA's
    """
    # get codon information into three lists
    gene_seq_str = gene_seq.exon_seq
    ref_codon, codon_pos, pos_in_codon = it.izip(*[cutils.pos_to_codon(gene_seq_str, p)
                                                   for p in coding_pos])
    ref_codon, codon_pos, pos_in_codon = list(ref_codon), list(codon_pos), list(pos_in_codon)

    # construct codons for mutations
    mut_codon = [list(x) for x in ref_codon]
    for i in range(len(mut_codon)):
        pc = pos_in_codon[i]
        mut_codon[i][pc] = somatic_base[i]
    mut_codon = [''.join(x) for x in mut_codon]

    # output resulting info
    aa_info = {'Reference Codon': ref_codon,
               'Somatic Codon': mut_codon,
               'Codon Pos': codon_pos,
               'Reference AA': [(codon_table[r] if len(r)==3 else None)
                                for r in ref_codon],
               'Somatic AA': [(codon_table[s] if len(s)==3 else None)
                              for s in mut_codon]}

    return aa_info
