cimport cython
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
import numpy as np
cimport numpy as np
from ..python import utils

# define data types for kde function
DTYPE_INT = np.int
# define compile time data types
ctypedef np.int_t DTYPE_INT_t


cdef extern from "permutation.hpp":
    # import functions from the C++ header permutation.hpp
    map[string, double] calc_position_statistics(map[int, int] pos_ctr,
                                                 double min_frac,
                                                 int min_recurrent,
                                                 int is_obs)
    map[string, double] calc_effect_statistics(map[int, int] pos_ctr,
                                               double min_frac,
                                               int min_recurrent,
                                               int is_obs)


@cython.cdivision(True)
def pos_to_codon(gene_seq, int pos):
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
    cdef int codon_pos, codon_start, pos_in_codon, seq_len = gene_seq.bed.cds_len
    if pos < seq_len:
        # valid mutation in coding region
        codon_pos = pos // 3
        codon_start = codon_pos * 3
        pos_in_codon = pos % 3
        ref = gene_seq.exon_seq[pos]
        return gene_seq.exon_seq[codon_start:codon_start+3], codon_pos, pos_in_codon, ref
    else:
        # by assumption, "positions" of splice sites are greater than the
        # length of the coding region to distinguish splice site mutations
        # from coding region mutations. To indicate the mutation is at a
        # splice site, I return None for positions.
        ss_pos = gene_seq.bed.pos2ss[pos]
        if ss_pos[0] == "5'":
            ref = gene_seq.five_prime_seq[ss_pos[1]][ss_pos[2]]
        else:
            ref = gene_seq.three_prime_seq[ss_pos[1]][ss_pos[2]]

        return 'Splice_Site', None, None, ref


def calc_pos_info(aa_mut_pos,
                  germ_aa,
                  somatic_aa,
                  int pseudo_count=0,
                  double min_frac=0.0,
                  int min_recur=2,
                  int is_obs=1):
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
    # pos_info = calc_position_statistics(pos_ctr, min_frac, min_recur, is_obs)
    pos_info = calc_position_statistics(pos_ctr, min_frac, min_recur, is_obs)
    num_recur = <int> pos_info["recurrent"]
    frac_pos_ent = pos_info["entropy_fraction"]
    delta_pos_ent = pos_info["delta_entropy"]
    return num_recur, frac_pos_ent, delta_pos_ent


def calc_effect_info(aa_mut_pos,
                     germ_aa,
                     somatic_aa,
                     int pseudo_count=0,
                     double min_frac=0.0,
                     int min_recur=2,
                     int is_obs=1):
    cdef:
        map[int, int] pos_ctr
        map[string, double] pos_info
        int num_recur = 0
        double frac_pos_ent = 0.0
        int i, num_pos
        DTYPE_INT_t[::1] pos_array
        int DUMMY_INT = 9999999  # dummy pos if prior used
        int INACTIVATING_INT = -1  # pos for inactivating mutations
    tmp_pos_list = []
    num_pos = len(aa_mut_pos)
    for i in range(num_pos):
        pos = aa_mut_pos[i]
        # make sure mutation is missense
        if germ_aa[i] and somatic_aa[i] and germ_aa[i] != '*' and \
           somatic_aa[i] != '*' and germ_aa[i] != somatic_aa[i] and \
           pos != 0:
            # should have a position, but if not skip it
            if pos is not None:
                if pos_ctr.count(pos) == 0:
                    pos_ctr[pos] = 0
                pos_ctr[pos] += 1
                tmp_pos_list.append(pos)
        elif ((germ_aa[i] == '*' or somatic_aa[i] == '*' or pos==0) and \
              (germ_aa[i] != somatic_aa[i])) or \
             (germ_aa[i] == 'Splice_Site' or somatic_aa[i] == 'Splice_Site'):
            # case for inactivating mutations
            if pos_ctr.count(INACTIVATING_INT) == 0:
                pos_ctr[INACTIVATING_INT] = 0
            pos_ctr[INACTIVATING_INT] += 1
            tmp_pos_list.append(INACTIVATING_INT)

    # add pseudo-counts if specified
    if pseudo_count:
        pos_ctr[DUMMY_INT] = pseudo_count

    # get entropy of effect statistics
    effect_info = calc_effect_statistics(pos_ctr, min_frac, min_recur, is_obs)
    num_recur = <int> effect_info["recurrent_sum"]
    frac_effect_ent = effect_info["entropy_fraction"]
    num_inactivating = effect_info["inactivating_sum"]
    return frac_effect_ent, num_recur, num_inactivating


def calc_deleterious_info(germ_aa, somatic_aa, codon_pos):
    cdef:
        int i, num_mutations = 0, num_deleterious = 0

    num_mutations = len(somatic_aa)
    if len(germ_aa) != num_mutations:
        raise ValueError('There should be equal number of germline and somatic bases')

    for i in range(num_mutations):
        if germ_aa[i] and somatic_aa[i] and \
           ((germ_aa[i] == '*' or somatic_aa[i] == '*' or codon_pos[i]==0) and \
            germ_aa[i] != somatic_aa[i]) or \
           somatic_aa[i] == 'Splice_Site':
            num_deleterious += 1

    return num_deleterious


def calc_non_silent_info(germ_aa, somatic_aa, codon_pos):
    cdef:
        int i, num_mutations = 0
        int num_non_silent = 0, num_silent = 0, num_nonsense = 0
        int num_loststop = 0, num_splice_site = 0, num_missense = 0
        int num_loststart = 0

    num_mutations = len(somatic_aa)
    if len(germ_aa) != num_mutations:
        raise ValueError('There should be equal number of germline and somatic bases')

    for i in range(num_mutations):
        if (germ_aa[i] and somatic_aa[i]) or somatic_aa[i] == 'Splice_Site' or germ_aa[i] == 'Splice_Site':
            # count nonsense
            if (somatic_aa[i] != germ_aa[i]) and somatic_aa[i] == '*':
                num_nonsense += 1
                num_non_silent += 1
            # count lost stop
            elif (somatic_aa[i] != germ_aa[i]) and germ_aa[i] == '*':
                num_loststop += 1
                num_non_silent += 1
            elif somatic_aa[i] == 'Splice_Site' or germ_aa[i] == 'Splice_Site':
                num_splice_site += 1
                num_non_silent += 1
            elif somatic_aa[i] != germ_aa[i]:
                if codon_pos[i] == 0:
                    num_loststart += 1
                else:
                    num_missense += 1
                num_non_silent += 1
            else:
                num_silent += 1

    return [num_non_silent, num_silent, num_nonsense,
            num_loststop, num_splice_site, num_loststart, num_missense]


def get_variant_classification(germ_aa_list, somatic_aa_list, codon_pos):
    """Return the proper variant classification for a substiution mutation.

    Parameters
    ----------
    germ_aa : list of str
        reference amino acid residue
    somatic_aa : list of str
        new mutated residues
    codon_pos : list of int
        codon position in sequence

    Returns
    -------
    var_class : list of str
        list of strings classifying variant type
    """
    cdef:
        string na = ''
        vector[string] var_class
        vector[string] germ_aa = [(g if g else na) for g in germ_aa_list]
        vector[string] somatic_aa = [(s if s else na) for s in somatic_aa_list]
        int num_muts = len(somatic_aa)
        string missense = 'Missense_Mutation'
        string nonsense = 'Nonsense_Mutation'
        string loststop = 'Nonstop_Mutation'
        string splice_site = 'Splice_Site'
        string silent = 'Silent'
        string lost_start = 'Translation_Start_Site'
        string stop_codon = '*'

    for i in range(num_muts):
        if (germ_aa[i].length() and somatic_aa[i].length()) or somatic_aa[i] == splice_site or germ_aa[i] == splice_site:
            # count nonsense
            if (somatic_aa[i] != germ_aa[i]) and somatic_aa[i] == stop_codon:
                var_class.push_back(nonsense)
            # count lost stop
            elif (somatic_aa[i] != germ_aa[i]) and germ_aa[i] == stop_codon:
                var_class.push_back(loststop)
            elif somatic_aa[i] == splice_site or germ_aa[i] == splice_site:
                var_class.push_back(splice_site)
            elif somatic_aa[i] != germ_aa[i]:
                if codon_pos[i] == 0:
                    var_class.push_back(lost_start)
                else:
                    var_class.push_back(missense)
            else:
                var_class.push_back(silent)
        else:
            var_class.push_back(na)

    return var_class


def calc_summary_info(germ_aa, somatic_aa, codon_pos,
                      min_frac=0.0,
                      min_recur=2):
    """Returns information from both missense position metrics and number of
    mutation types.

    Mostly a wrapper around calc_non_silent_info and calc_pos_info.

    Parameters
    ----------
    germ_aa : list
        list of strings indicating reference amino acid (single letter)
    somatic_aa : list
        list of strings indicating mutated amino acid
    codon_pos : list
        contains integer position of codon
    min_frac : float
        fraction of total mutations to be recurrent position
    min_recur : int
        minimum number of missense at same position to be defined as recurrent

    Returns
    -------
    summary information
    """
    mut_type_info = calc_non_silent_info(germ_aa, somatic_aa, codon_pos)
    num_recur, pos_ent, delta_ent = calc_pos_info(codon_pos, germ_aa,
                                                  somatic_aa,
                                                  min_frac=min_frac,
                                                  min_recur=min_recur
                                                  #is_obs=0
                                                  )
    return mut_type_info + [num_recur, pos_ent]
