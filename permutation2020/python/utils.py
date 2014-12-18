# normal imports
from bed_line import BedLine
import numpy as np
import scipy.stats as stats
import pandas as pd
import csv
import itertools as it
from functools import wraps
import logging

logger = logging.getLogger(__name__)  # module logger

# small epsilon value to prevent issues with machine decimal precision
epsilon = 0.0001

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

def log_error_decorator(f):
    """Writes exception to log file if occured in decorated function.

    This decorator wrapper is needed for multiprocess logging since otherwise
    the python multiprocessing module will obscure the actual line of the error.
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        try:
            result = f(*args, **kwds)
            return result
        except KeyboardInterrupt:
            logger.info('Ctrl-C stopped a process.')
        except Exception, e:
            logger.exception(e)
            raise
    return wrapper


def keyboard_exit_wrapper(func):
    def wrap(self, timeout=None):
        # Note: the timeout of 1 googol seconds introduces a rather subtle
        # bug for Python scripts intended to run many times the age of the universe.
        return func(self, timeout=timeout if timeout is not None else 1e100)
    return wrap


def filter_list(mylist, bad_ixs):
    """Removes indices from a list.

    All elements in bad_ixs will be removed from the list.

    Parameters
    ----------
    mylist : list
        list to filter out specific indices
    bad_ixs : list of ints
        indices to remove from list

    Returns
    -------
    mylist : list
        list with elements filtered out
    """
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
    """Check if valid single letter base.

    Parameters
    ----------
    nuc : str
        String to check if valid single nucleotide base

    Returns
    -------
    is_valid : bool
        flag indicating valid nucleotide base
    """
    valid_nucs = ['A', 'C', 'T', 'G', 'N']
    is_valid = nuc in valid_nucs
    return is_valid


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


def read_bed(file_path, filtered_genes):
    """Reads BED file and populates a dictionary separating genes
    by chromosome.

    Parameters
    ----------
    file_path : str
        path to BED file
    filtered_genes: list
        list of gene names to not use

    Returns
    -------
    bed_dict: dict
        dictionary mapping chromosome keys to a list of BED lines
    """
    # read in entire bed file into a dict with keys as chromsomes
    bed_dict = {}
    for bed_row in bed_generator(file_path):
        if bed_row.gene_name not in filtered_genes:
            bed_dict.setdefault(bed_row.chrom, [])
            bed_dict[bed_row.chrom].append(bed_row)
    return bed_dict


def _fix_mutation_df(mutation_df):
    """Drops invalid mutations and corrects for 1-based coordinates.

    TODO: Be smarter about what coordinate system is put in the provided
    mutations.

    Parameters
    ----------
    mutation_df : pd.DataFrame
        user provided mutations

    Returns
    -------
    mutation_df : pd.DataFrame
        mutations filtered for being valid and correct mutation type. Also
        converted 1-base coordinates to 0-based.

    """
    # only keep allowed mutation types
    orig_len = len(mutation_df)  # number of mutations before filtering
    allowed_types = ['Missense_Mutation', 'Silent', 'Nonsense_Mutation',
                     'Splice_Site', 'Nonstop_Mutation']
    mutation_df = mutation_df[mutation_df.Variant_Classification.isin(allowed_types)]  # only keep SNV
    type_len = len(mutation_df)  # number of mutations after filtering based on mut type

    # log the number of dropped mutations
    log_msg = ('Dropped {num_dropped} mutations after only keeping '
               '{mut_types}'.format(num_dropped=orig_len-type_len,
                                    mut_types=', '.join(allowed_types)))
    logger.info(log_msg)

    # check if mutations are valid SNVs
    valid_nuc_flag = (mutation_df['Reference_Allele'].apply(is_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(is_valid_nuc))
    mutation_df = mutation_df[valid_nuc_flag]  # filter bad lines
    mutation_df = mutation_df[mutation_df['Tumor_Allele'].apply(lambda x: len(x)==1)]
    mutation_df = mutation_df[mutation_df['Reference_Allele'].apply(lambda x: len(x)==1)]
    valid_len = len(mutation_df)

    # log the number of dropped mutations
    log_msg = ('Dropped {num_dropped} mutations after only keeping '
               'valid SNVs'.format(num_dropped=type_len-valid_len))
    logger.info(log_msg)

    # correct for 1-based coordinates
    mutation_df['Start_Position'] = mutation_df['Start_Position'] - 1
    return mutation_df


def _get_high_tsg_score(mutation_df, tsg_score_thresh):
    # find genes above a tsg score threshold
    mutation_df['indicator'] = 1
    table = pd.pivot_table(mutation_df,
                           values='indicator',
                           cols='Variant_Classification',
                           rows='Gene',
                           aggfunc=np.sum)
    mut_type_frac = table.div(table.sum(axis=1).astype(float), axis=0).fillna(0.0)
    for c in ['Nonsense_Mutation', 'Frame_Shift_Indel', 'Splice_Site', 'Nonstop_Mutation']:
        if c not in mut_type_frac.columns:
            mut_type_frac[c] = 0.0  # make sure columns are defined
    tsg_score = mut_type_frac['Nonsense_Mutation'] + mut_type_frac['Frame_Shift_Indel'] + \
                mut_type_frac['Splice_Site'] + mut_type_frac['Nonstop_Mutation']
    non_tested_genes = set(tsg_score[tsg_score>=tsg_score_thresh].index.tolist())

    # log the number of non tested genes
    log_msg = ('{0} genes will not be tested due to high TSG '
               'score'.format(len(non_tested_genes)))
    logger.info(log_msg)

    return non_tested_genes


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


def fishers_method(pvals):
    pvals = np.asarray(pvals)
    degrees_of_freedom = 2 * pvals.size
    chisq_stat = np.sum(-2*np.log(pvals))
    fishers_pval = stats.chi2.sf(chisq_stat, degrees_of_freedom)
    return fishers_pval


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


