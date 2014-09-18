cimport cython
from libcpp.map cimport map
from libcpp.string cimport string
import numpy as np
cimport numpy as np


# define data types for kde function
DTYPE_INT = np.int
# define compile time data types
ctypedef np.int_t DTYPE_INT_t


cdef extern from "permutation.hpp":
    # import functions from the C++ header permutation.hpp
    map[string, double] calc_position_statistics(map[int, int] pos_ctr, double min_frac)


@cython.cdivision(True)
def pos_to_codon(seq, int pos):
    """Retrieves information about the codon a nucleotide position is in.

    Parameters
    ----------
    seq : str
        coding sequence
    pos : int
        0-based position of nucleotide in seq

    Returns
    -------
    seq : str
        actual codon sequence
    codon_pos : int
        0-based position of codon (e.g. 3 is the 4th codon)
    pos_in_codon : int
        0-based position within a codon (e.g. 1 is the second
        position out of three)
    """
    cdef int codon_pos, codon_start, pos_in_codon, seq_len = len(seq)
    if pos < seq_len:
        # valid mutation in coding region
        codon_pos = pos / 3
        codon_start = codon_pos * 3
        pos_in_codon = pos % 3
        return seq[codon_start:codon_start+3], codon_pos, pos_in_codon
    else:
        # by assumption, "positions" of splice sites are greater than the
        # length of the coding region to distinguish splice site mutations
        # from coding region mutations. To indicate the mutation is at a
        # splice site, I return None for positions.
        return 'Splice_Site', None, None


def calc_pos_info(aa_mut_pos,
                  germ_aa,
                  somatic_aa,
                  int pseudo_count=0,
                  double min_frac=0.0):
    cdef:
        map[int, int] pos_ctr
        map[string, double] pos_info
        int num_recur = 0
        double frac_pos_ent = 0.0
        double delta_pos_ent = 0.0
        int i, num_pos
        DTYPE_INT_t[::1] pos_array
        cdef int DUMMY_INT = 9999999  # dummy pos if prior used
    tmp_pos_list = []
    num_pos = len(aa_mut_pos)
    for i in range(num_pos):
        pos = aa_mut_pos[i]
        # make sure mutation is missense
        if germ_aa[i] and somatic_aa[i] and germ_aa[i] != '*' and \
           somatic_aa[i] != '*' and germ_aa[i] != somatic_aa[i]:
            # should have a position, but if not skip it
            if pos is not None:
                if pos_ctr.count(pos) == 0:
                    pos_ctr[pos] = 0
                pos_ctr[pos] += 1
                tmp_pos_list.append(pos)

    # add pseudo-counts if specified
    if pseudo_count:
        pos_ctr[DUMMY_INT] = pseudo_count

    # get position statistics
    pos_info = calc_position_statistics(pos_ctr, min_frac)
    num_recur = <int> pos_info["recurrent"]
    frac_pos_ent = pos_info["entropy_fraction"]
    delta_pos_ent = pos_info["delta_entropy"]
    return num_recur, frac_pos_ent, delta_pos_ent


def calc_deleterious_info(germ_aa, somatic_aa):
    cdef:
        int i, num_mutations = 0, num_deleterious = 0

    num_mutations = len(somatic_aa)
    if len(germ_aa) != num_mutations:
        raise ValueError('There should be equal number of germline and somatic bases')

    for i in range(num_mutations):
        if germ_aa[i] and somatic_aa[i] and \
           ((germ_aa[i] == '*' or somatic_aa[i] == '*') and \
            germ_aa[i] != somatic_aa[i]) or \
           somatic_aa[i] == 'Splice_Site':
            num_deleterious += 1

    return num_deleterious


def calc_non_silent_info(germ_aa, somatic_aa):
    cdef:
        int i, num_mutations = 0, num_non_silent = 0, num_silent = 0

    num_mutations = len(somatic_aa)
    if len(germ_aa) != num_mutations:
        raise ValueError('There should be equal number of germline and somatic bases')

    for i in range(num_mutations):
        if germ_aa[i] and somatic_aa[i] and \
           somatic_aa[i] == 'Splice_Site' or \
           somatic_aa[i] != germ_aa[i]:
            num_non_silent += 1
        else:
            num_silent += 1

    return [num_non_silent, num_silent]
