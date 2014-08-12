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
from amino_acid import AminoAcid


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
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*',
               'Splice_Site': 'Splice_Site'}

# global dictionary specifying base pairing
base_pairing = {'A': 'T',
                'T': 'A',
                'a': 't',
                't': 'a',
                'C': 'G',
                'G': 'C',
                'c': 'g',
                'g': 'c',
                '-': '-',  # some people denote indels with '-'
                'n': 'n',
                'N': 'N'}

def filter_list(mylist, bad_ixs):
    # indices need to be in reverse order for filtering
    # to prevent .pop() from yielding eroneous results
    bad_ixs = sorted(bad_ixs, reverse=True)
    for i in bad_ixs:
        mylist.pop(i)
    return mylist


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
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base_pairing[s] for s in rev_seq])
    return rev_comp_seq


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


def get_chasm_context(tri_nuc):
    """Returns the mutation context acording to CHASM.

    For more information about CHASM's mutation context, look
    at http://wiki.chasmsoftware.org/index.php/CHASM_Overview.
    Essentially CHASM uses a few specified di-nucleotide contexts
    followed by single nucleotide context.

    Parameters
    ----------
    tri_nuc : str
        three nucleotide string with mutated base in the middle.

    Returns
    -------
    chasm context : str
        a string representing the context used in CHASM
    """
    # check if string is correct length
    if len(tri_nuc) != 3:
        raise ValueError('Chasm context requires a three nucleotide string '
                         '(Provided: "{0}")'.format(tri_nuc))

    # try dinuc context if found
    if tri_nuc[1:] == 'CG':
        return 'C*pG'
    elif tri_nuc[:2] == 'CG':
        return 'CpG*'
    elif tri_nuc[:2] == 'TC':
        return 'TpC*'
    elif tri_nuc[1:] == 'GA':
        return 'G*pA'
    else:
        # just return single nuc context
        return tri_nuc[1]


def get_context(chr, pos_list, strand, fa, context_type):
    nuc_contexts = []
    if context_type in [1, 2]:
        # case where context matters
        index_context = int(context_type) - 1  # subtract 1 since python is zero-based index
        for pos in pos_list:
            nucs = fa.fetch(reference=chr,
                            start=pos-index_context,
                            end=pos+1).upper()
            if strand == '-':
                nucs = rev_comp(nucs)
            if 'N' not in nucs:
                nuc_contexts.append(nucs)
            else:
                nuc_contexts.append(None)
    elif context_type in [1.5, 3]:
        # use the nucleotide context from chasm if nuc
        # context is 1.5 otherwise always use a three
        # nucleotide context
        for pos in pos_list:
            nucs = fa.fetch(reference=chr,
                            start=pos-1,
                            end=pos+2).upper()
            if strand == '-':
                nucs = rev_comp(nucs)
            if context_type == 1.5 and nucs:
                nucs = get_chasm_context(nucs)

            if 'N' not in nucs:
                nuc_contexts.append(nucs)
            else:
                nuc_contexts.append(None)
    else:
        nuc_contexts = ['None'] * len(pos_list)

    return nuc_contexts


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
    mut_codon = [(list(x) if x != 'Splice_Site' else []) for x in ref_codon]
    for i in range(len(mut_codon)):
        # splice site mutations are not in a codon, so skip such mutations to
        # prevent an error
        if pos_in_codon[i] is not None:
            pc = pos_in_codon[i]
            mut_codon[i][pc] = somatic_base[i]
    mut_codon = [(''.join(x) if x else 'Splice_Site') for x in mut_codon]

    # output resulting info
    aa_info = {'Reference Codon': ref_codon,
               'Somatic Codon': mut_codon,
               'Codon Pos': codon_pos,
               'Reference AA': [(codon_table[r] if len(r)==3 else None)
                                for r in ref_codon],
               'Somatic AA': [(codon_table[s] if len(s)==3 else None)
                              for s in mut_codon]}

    return aa_info


def get_unmapped_aa_mut_info(mut_info, genome_fa, strand, chr, context_type):

    # get information on the nucleotide context
    mycontexts = get_context(chr, mut_info['Start_Position'],
                             strand, genome_fa, context_type)

    # get information about the effect of the protein change
    if len(mut_info) > 0:
        not_splice_site = mut_info['Variant_Classification'].map(lambda x: x!='Splice_Site')
        prot_change = [AminoAcid(p) for p in mut_info[not_splice_site]['Protein_Change']]
    else:
        prot_change = []
    codon_pos, germ_aa, somatic_aa = [], [], []
    LARGE_NUMBER = 100000  # sufficiently large number to prevent accidental overlap of codon positions
    tmp_index = 0
    bad_mut_ix, good_mut_ix = [], []
    for i in range(len(mut_info)):
        if not mycontexts[i]:
            bad_mut_ix.append(i)  # remove invalid/missing mutation
            codon_pos.append(None)
            germ_aa.append(None)
            somatic_aa.append(None)
            if not_splice_site.iloc[i]:
                tmp_index += 1
        elif not_splice_site.iloc[i]:
            if prot_change and (not prot_change[tmp_index].is_valid or \
                                prot_change[tmp_index].is_missing_info):
                bad_mut_ix.append(i)  # remove invalid/missing mutation
                codon_pos.append(None)
                germ_aa.append(None)
                somatic_aa.append(None)
                tmp_index += 1
            else:
                good_mut_ix.append(i)
                codon_pos.append(LARGE_NUMBER + prot_change[tmp_index].pos)
                germ_aa.append(prot_change[tmp_index].initial)
                somatic_aa.append(prot_change[tmp_index].mutated)
                tmp_index += 1
        else:
            good_mut_ix.append(i)
            codon_pos.append('Splice_Site')
            germ_aa.append('Splice_Site')
            somatic_aa.append('Splice_Site')

    # remove bad mutations from results
    mycontexts = filter_list(mycontexts, bad_mut_ix)
    germ_aa = filter_list(germ_aa, bad_mut_ix)
    somatic_aa = filter_list(somatic_aa, bad_mut_ix)
    codon_pos = filter_list(codon_pos, bad_mut_ix)

    # information about the effect of mutations that could not be mapped
    # to the reference isoform of a gene.
    tumor_allele = mut_info.iloc[np.array(good_mut_ix, dtype=int)]['Tumor_Allele'].tolist()
    aa_info = {'Context': mycontexts,
               'Codon Pos': codon_pos,
               'Reference AA': germ_aa,
               'Somatic AA': somatic_aa,
               'Tumor_Allele': tumor_allele}

    return aa_info
