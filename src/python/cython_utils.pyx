from libcpp.map cimport map

cdef extern from "permutation.hpp":
    int recurrent_sum(map[int, int] pos_ctr)
    double position_entropy(map[int, int] pos_ctr)


def pos_to_codon(seq, int pos):
    cdef int codon_pos, codon_start, pos_in_codon
    codon_pos = pos / 3
    codon_start = codon_pos * 3
    pos_in_codon = pos % 3
    return seq[codon_start:codon_start+3], codon_pos, pos_in_codon


def calc_pos_info(aa_mut_pos):
    cdef:
        map[int, int] pos_ctr
        int num_recur = 0
        double pos_ent = 0.0
    for pos in aa_mut_pos:
        if pos is not None:
            if pos_ctr.count(pos) == 0:
                pos_ctr[pos] = 0
            pos_ctr[pos] += 1
    num_recur = recurrent_sum(pos_ctr)
    pos_ent = position_entropy(pos_ctr)
    return num_recur, pos_ent
