import numpy as np
import utils
from ..cython import cutils
import mutation_context as mc


def deleterious_permutation(context_counts,
                            context_to_mut,
                            seq_context,
                            gene_seq,
                            num_permutations=10000,
                            pseudo_count=0):
    """Performs null-permutations for deleterious mutation statistics
    in a single gene.

    Parameters
    ----------
    context_counts : pd.Series
        number of mutations for each context
    context_to_mut : dict
        dictionary mapping nucleotide context to a list of observed
        somatic base changes.
    seq_context : SequenceContext
        Sequence context for the entire gene sequence (regardless
        of where mutations occur). The nucleotide contexts are
        identified at positions along the gene.
    gene_seq : GeneSequence
        Sequence of gene of interest
    num_permutations : int, default: 10000
        number of permutations to create for null
    pseudo_count : int, default: 0
        Pseudo-count for number of deleterious mutations for each
        permutation of the null distribution. Increasing pseudo_count
        makes the statistical test more stringent.

    Returns
    -------
    del_count_list : list
        list of deleterious mutation counts under the null
    """
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack(pos_array for base, pos_array in tmp_contxt_pos)

    # determine result of random positions
    del_count_list = []
    for row in tmp_mut_pos:
        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # calc deleterious mutation info
        tmp_del_count = cutils.calc_deleterious_info(tmp_mut_info['Reference AA'],
                                                     tmp_mut_info['Somatic AA'],
                                                     tmp_mut_info['Codon Pos'])
        del_count_list.append(tmp_del_count + pseudo_count)
    return del_count_list


def position_permutation(context_counts,
                         context_to_mut,
                         seq_context,
                         gene_seq,
                         num_permutations=10000,
                         pseudo_count=0):
    """Performs null-permutations for position-based mutation statistics
    in a single gene.

    Parameters
    ----------
    context_counts : pd.Series
        number of mutations for each context
    context_to_mut : dict
        dictionary mapping nucleotide context to a list of observed
        somatic base changes.
    seq_context : SequenceContext
        Sequence context for the entire gene sequence (regardless
        of where mutations occur). The nucleotide contexts are
        identified at positions along the gene.
    gene_seq : GeneSequence
        Sequence of gene of interest
    num_permutations : int, default: 10000
        number of permutations to create for null
    pseudo_count : int, default: 0
        Pseudo-count for number of recurrent missense mutations for each
        permutation for the null distribution. Increasing pseudo_count
        makes the statistical test more stringent.

    Returns
    -------
    num_recur_list : list
        list of recurrent mutation counts under the null
    entropy_list : list
        list of position entropy values under the null
    """
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack(pos_array for base, pos_array in tmp_contxt_pos)

    # calculate position-based statistics as a result of random positions
    #num_recur_list, entropy_list, kde_entropy_list, bw_list = [], [], [], []
    num_recur_list, entropy_list, delta_entropy_list = [], [], []
    for row in tmp_mut_pos:
        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # calculate position info
        tmp_recur_ct, tmp_entropy, tmp_delta_entropy = cutils.calc_pos_info(tmp_mut_info['Codon Pos'],
                                                                            tmp_mut_info['Reference AA'],
                                                                            tmp_mut_info['Somatic AA'],
                                                                            pseudo_count=pseudo_count,
                                                                            is_obs=0)
        num_recur_list.append(tmp_recur_ct)
        entropy_list.append(tmp_entropy)
        delta_entropy_list.append(tmp_delta_entropy)

    return num_recur_list, entropy_list, delta_entropy_list


def effect_permutation(context_counts,
                       context_to_mut,
                       seq_context,
                       gene_seq,
                       num_permutations=10000,
                       pseudo_count=0):
    """Performs null-permutations for effect-based mutation statistics
    in a single gene.

    Parameters
    ----------
    context_counts : pd.Series
        number of mutations for each context
    context_to_mut : dict
        dictionary mapping nucleotide context to a list of observed
        somatic base changes.
    seq_context : SequenceContext
        Sequence context for the entire gene sequence (regardless
        of where mutations occur). The nucleotide contexts are
        identified at positions along the gene.
    gene_seq : GeneSequence
        Sequence of gene of interest
    num_permutations : int, default: 10000
        number of permutations to create for null
    pseudo_count : int, default: 0
        Pseudo-count for number of recurrent missense mutations for each
        permutation for the null distribution. Increasing pseudo_count
        makes the statistical test more stringent.

    Returns
    -------
    effect_entropy_list : list
        list of entropy of effect values under the null
    recur_list : list
        number of recurrent missense mutations
    inactivating_list : list
        number of inactivating mutations
    """
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack(pos_array for base, pos_array in tmp_contxt_pos)

    # calculate position-based statistics as a result of random positions
    effect_entropy_list, recur_list, inactivating_list = [], [], []
    for row in tmp_mut_pos:
        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # calculate position info
        tmp_entropy, tmp_recur, tmp_inactivating = cutils.calc_effect_info(tmp_mut_info['Codon Pos'],
                                                                           tmp_mut_info['Reference AA'],
                                                                           tmp_mut_info['Somatic AA'],
                                                                           pseudo_count=pseudo_count,
                                                                           is_obs=0)
        effect_entropy_list.append(tmp_entropy)
        recur_list.append(tmp_recur)
        inactivating_list.append(tmp_inactivating)

    return effect_entropy_list, recur_list, inactivating_list


def non_silent_ratio_permutation(context_counts,
                                 context_to_mut,
                                 seq_context,
                                 gene_seq,
                                 num_permutations=10000):
    """Performs null-permutations for non-silent ratio across all genes.

    Parameters
    ----------
    context_counts : pd.Series
        number of mutations for each context
    context_to_mut : dict
        dictionary mapping nucleotide context to a list of observed
        somatic base changes.
    seq_context : SequenceContext
        Sequence context for the entire gene sequence (regardless
        of where mutations occur). The nucleotide contexts are
        identified at positions along the gene.
    gene_seq : GeneSequence
        Sequence of gene of interest
    num_permutations : int, default: 10000
        number of permutations to create for null

    Returns
    -------
    non_silent_count_list : list of tuples
        list of non-silent and silent mutation counts under the null
    """
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack(pos_array for base, pos_array in tmp_contxt_pos)

    # determine result of random positions
    non_silent_count_list = []
    for row in tmp_mut_pos:
        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # calc deleterious mutation info
        tmp_non_silent = cutils.calc_non_silent_info(tmp_mut_info['Reference AA'],
                                                     tmp_mut_info['Somatic AA'],
                                                     tmp_mut_info['Codon Pos'])
        non_silent_count_list.append(tmp_non_silent)
    return non_silent_count_list
