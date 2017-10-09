# normal imports
from prob2020.python.bed_line import BedLine
import numpy as np
import pandas as pd
import csv
from collections import OrderedDict
from functools import wraps
import warnings

# logging import
import logging
import datetime
import os
import sys

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

##############################
# Define groups of mutation consequences
# based on the variant classification
# column names
##############################
# variant classifications for SNVs
variant_snv = ['Missense_Mutation', 'Silent', 'Nonsense_Mutation',
               'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site']

# variant classifications that are inactivating
variant_inactivating = ['Splice_Site', 'Nonsense_Mutation', 'Translation_Start_Site',
                        'Nonstop_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
                        'Frame_Shift_Indel']

# coding variant classifications that are not necessarily inactivating
variant_non_inactivating = ['Missense_Mutation', 'Silent']

# indel/frameshift variant names
variant_frameshift = ['Frame_Shift_Indel', 'Frame_Shift_Ins', 'Frame_Shift_Del']
variant_in_frame_indel = ['In_Frame_Indel', 'In_Frame_Ins', 'In_Frame_Del']
variant_indel = variant_frameshift + variant_in_frame_indel

# all variants
all_variants = variant_snv + variant_indel

def start_logging(log_file='', log_level='INFO', verbose=False):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO

    # ignore warnings if not in debug
    if log_level.upper() != 'DEBUG':
        warnings.filterwarnings('ignore')

    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

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
        except Exception as e:
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


def read_bed(file_path, restricted_genes=None):
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
    bed_dict = OrderedDict()
    for bed_row in bed_generator(file_path):
        is_restrict_flag = restricted_genes is None or bed_row.gene_name in restricted_genes
        if is_restrict_flag:
            bed_dict.setdefault(bed_row.chrom, [])
            bed_dict[bed_row.chrom].append(bed_row)
    sort_chroms = sorted(bed_dict.keys(), key=lambda x: len(bed_dict[x]), reverse=True)
    bed_dict = OrderedDict((chrom, bed_dict[chrom]) for chrom in sort_chroms)
    return bed_dict


def _fix_mutation_df(mutation_df, only_unique=False):
    """Drops invalid mutations and corrects for 1-based coordinates.

    TODO: Be smarter about what coordinate system is put in the provided
    mutations.

    Parameters
    ----------
    mutation_df : pd.DataFrame
        user provided mutations
    only_unique : bool
        flag indicating whether only unique mutations for each tumor sample
        should be kept. This avoids issues when the same mutation has
        duplicate reportings.

    Returns
    -------
    mutation_df : pd.DataFrame
        mutations filtered for being valid and correct mutation type. Also
        converted 1-base coordinates to 0-based.

    """
    # only keep allowed mutation types
    orig_len = len(mutation_df)  # number of mutations before filtering
    mutation_df = mutation_df[mutation_df.Variant_Classification.isin(variant_snv)]  # only keep SNV
    type_len = len(mutation_df)  # number of mutations after filtering based on mut type

    # log the number of dropped mutations
    log_msg = ('Dropped {num_dropped} mutations after only keeping '
               '{mut_types}. Indels are processed separately.'.format(num_dropped=orig_len-type_len,
                                                                      mut_types=', '.join(variant_snv)))
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

    # drop duplicate mutations
    if only_unique:
        dup_cols = ['Tumor_Sample', 'Chromosome', 'Start_Position',
                    'End_Position', 'Reference_Allele', 'Tumor_Allele']
        mutation_df = mutation_df.drop_duplicates(subset=dup_cols)

        # log results of de-duplication
        dedup_len = len(mutation_df)
        log_msg = ('Dropped {num_dropped} mutations when removing '
                   'duplicates'.format(num_dropped=valid_len-dedup_len))
        logger.info(log_msg)

    # add dummy Protein_Change or Tumor_Type columns if not provided
    # in file
    if 'Tumor_Type' not in mutation_df.columns:
        mutation_df['Tumor_Type'] = ''
    if 'Protein_Change' not in mutation_df.columns:
        mutation_df['Protein_Change'] = ''

    # correct for 1-based coordinates
    mutation_df['Start_Position'] = mutation_df['Start_Position'].astype(int) - 1
    return mutation_df


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


def calc_windowed_sum(aa_mut_pos,
                      germ_aa,
                      somatic_aa,
                      window=[3]):
    """Calculate the sum of mutations within a window around a particular mutated
    codon.

    Parameters
    ----------
    aa_mut_pos : list
        list of mutated amino acid positions
    germ_aa : list
        Reference amino acid
    somatic_aa : list
        Somatic amino acid (if missense)
    window : list
        List of windows to calculate for

    Returns
    -------
    pos_ctr : dict
        dictionary of mutated positions (key) with associated counts (value)
    pos_sum : dict of dict
        Window size as first key points to dictionary of mutated positions (key)
        with associated mutation count within the window size (value)
    """
    pos_ctr, pos_sum = {}, {w: {} for w in window}
    num_pos = len(aa_mut_pos)
    # figure out the missense mutations
    for i in range(num_pos):
        pos = aa_mut_pos[i]
        # make sure mutation is missense
        if germ_aa[i] and somatic_aa[i] and germ_aa[i] != '*' and \
           somatic_aa[i] != '*' and germ_aa[i] != somatic_aa[i]:
            # should have a position, but if not skip it
            if pos is not None:
                pos_ctr.setdefault(pos, 0)
                pos_ctr[pos] += 1

    # calculate windowed sum
    pos_list = sorted(pos_ctr.keys())
    max_window = max(window)
    for ix, pos in enumerate(pos_list):
        tmp_sum = {w: 0 for w in window}
        # go through the same and lower positions
        for k in reversed(range(ix+1)):
            pos2 = pos_list[k]
            if pos2 < pos-max_window: break
            for w in window:
                if pos-w <= pos2:
                    tmp_sum[w] += pos_ctr[pos2]
        # go though the higher positions
        for l in range(ix+1, len(pos_list)):
            pos2 = pos_list[l]
            if pos2 > pos+max_window: break
            for w in window:
                if pos2 <= pos+w:
                    tmp_sum[w] += pos_ctr[pos2]
        # iterate through all other positions
        #for pos2 in pos_list:
            #for w in window:
                #if pos-w <= pos2 <= pos+w:
                    #tmp_sum[w] += pos_ctr[pos2]
        # update windowed counts
        for w in window:
            pos_sum[w][pos] = tmp_sum[w]

    return pos_ctr, pos_sum
