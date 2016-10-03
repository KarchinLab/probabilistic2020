#from ..cython import cutils
import numpy as np
import os
import prob2020.python.mymath as mymath
import sys

# import pickle module
try:
    import cPickle as pickle
except:
    if sys.version_info < (3,):
        print('Falling back to regular pickle module')
    import pickle as pickle

def retrieve_scores(gname, sdir,
                    codon_pos, germ_aa, somatic_aa,
                    default_mga=5., default_vest=0,
                    no_file_flag=-1):
    """Retrieves scores from pickle files.

    Used by summary script.

    """
    # get variant types
    #var_class = cutils.get_variant_classification(germ_aa, somatic_aa, codon_pos)

    # get information about MGA entropy
    mga_path = os.path.join(sdir, gname+".mgaentropy.pickle")
    if os.path.exists(mga_path):
        if sys.version_info < (3,):
            # python 2.7 way
            with open(mga_path) as handle:
                mga_ent = pickle.load(handle)
        else:
            # python 3.X way
            with open(mga_path, 'rb') as handle:
                mga_ent = pickle.load(handle, encoding='latin-1')
    else:
        mga_ent = None
    missense_pos = [p for i, p in enumerate(codon_pos)
                    if (germ_aa[i]!=somatic_aa[i]) and
                       (germ_aa[i] not in ['-', '*', 'Splice_Site']) and
                       (somatic_aa[i] not in ['-', '*', 'Splice_Site'])]
    total_mga_ent = compute_mga_entropy_stat(mga_ent, missense_pos, sum, default_mga)

        #mga_ent_ixs = [codon_pos[i] for i in range(len(var_class))
                       #if var_class[i] == 'Missense_Mutation']
        #len_mga_ent = len(mga_ent)
        #mga_ent_scores = [mga_ent[ix] for ix in mga_ent_ixs if ix < len_mga_ent]
        #if mga_ent_scores:
            #total_mga_ent = sum(mga_ent_scores)
        #else:
            #total_mga_ent = default_mga
    #else:
        #total_mga_ent = no_file_flag

    # get information about VEST scores
    vest_path = os.path.join(sdir, gname+".vest.pickle")
    if os.path.exists(vest_path):
        if sys.version_info < (3,):
            # python 2.7 way
            with open(vest_path) as handle:
                vest_score = pickle.load(handle)
        else:
            # python 3.X way
            with open(vest_path, 'rb') as handle:
                vest_score = pickle.load(handle, encoding='latin-1')
    else:
        vest_score = None
    total_vest = compute_vest_stat(vest_score,
                                   germ_aa, somatic_aa, codon_pos,
                                   stat_func=sum, default_val=default_vest)
        #vest_scores = [vest_score.get(codon_pos[i]+1, {}).get(germ_aa[i], {}).get(somatic_aa[i], default_vest)
                        #for i in range(len(var_class))
                        #if var_class[i] == 'Missense_Mutation']
        #total_vest = sum(vest_scores)
    #else:
        #total_vest = no_file_flag
    return total_mga_ent, total_vest


def read_vest_pickle(gname, score_dir):
    """Read in VEST scores for given gene.

    Parameters
    ----------
    gname : str
        name of gene
    score_dir : str
        directory containing vest scores

    Returns
    -------
    gene_vest : dict or None
        dict containing vest scores for gene. Returns None if not found.
    """
    vest_path = os.path.join(score_dir, gname+".vest.pickle")
    if os.path.exists(vest_path):
        if sys.version_info < (3,):
            with open(vest_path) as handle:
                gene_vest = pickle.load(handle)
        else:
            with open(vest_path, 'rb') as handle:
                gene_vest = pickle.load(handle, encoding='latin-1')
        return gene_vest
    else:
        return None


def compute_vest_stat(vest_dict, ref_aa, somatic_aa, codon_pos,
                      stat_func=np.mean,
                      default_val=0.0):
    """Compute missense VEST score statistic.

    Note: non-missense mutations are intentially not filtered out and will take
    a default value of zero.

    Parameters
    ----------
    vest_dict : dict
        dictionary containing vest scores across the gene of interest
    ref_aa: list of str
        list of reference amino acids
    somatic_aa: list of str
        somatic mutation aa
    codon_pos : list of int
        position of codon in protein sequence
    stat_func : function, default=np.mean
        function that calculates a statistic
    default_val : float
        default value to return if there are no mutations

    Returns
    -------
    score_stat : float
        vest score statistic for provided mutation list
    """
    # return default value if VEST scores are missing
    if vest_dict is None:
        return default_val

    # fetch scores
    myscores = fetch_vest_scores(vest_dict, ref_aa, somatic_aa, codon_pos)

    # calculate mean score
    if myscores:
        score_stat = stat_func(myscores)
    else:
        score_stat = default_val

    return score_stat


def compute_mga_entropy_stat(mga_vec, codon_pos,
                             stat_func=np.mean,
                             default_val=0.0):
    """Compute MGA entropy conservation statistic

    Parameters
    ----------
    mga_vec : np.array
        numpy vector containing MGA Entropy conservation scores for residues
    codon_pos : list of int
        position of codon in protein sequence
    stat_func : function, default=np.mean
        function that calculates a statistic
    default_val : float
        default value to return if there are no mutations

    Returns
    -------
    score_stat : float
        MGA entropy score statistic for provided mutation list
    """
    # return default value if VEST scores are missing
    if mga_vec is None:
        return default_val

    # fetch scores
    myscores = fetch_mga_scores(mga_vec, codon_pos)

    # calculate mean score
    if myscores is not None and len(myscores):
        score_stat = stat_func(myscores)
    else:
        score_stat = default_val

    return score_stat


def fetch_vest_scores(vest_dict,
                      ref_aa, somatic_aa, codon_pos,
                      default_vest=0.0):
    """Get VEST scores from pre-computed scores in dictionary.

    Note: either all mutations should be missense or non-missense intended
    to have value equal to default.

    Parameters
    ----------
    vest_dict : dict
        dictionary containing vest scores across the gene of interest
    ref_aa: list of str
        list of reference amino acids
    somatic_aa: list of str
        somatic mutation aa
    codon_pos: list of int
        position of codon in protein sequence
    default_vest: float, default=0.0
        value to use if VEST score not available for a given mutation

    Returns
    -------
    vest_score_list: list
        score results for mutations
    """
    vest_score_list = []
    for i in range(len(somatic_aa)):
        # make sure position is valid
        if codon_pos[i] is not None:
            tmp_score = vest_dict.get(codon_pos[i]+1, {}).get(ref_aa[i], {}).get(somatic_aa[i], default_vest)
        else:
            tmp_score = 0.0
        vest_score_list.append(tmp_score)
    return vest_score_list


def fetch_mga_scores(mga_vec,
                     codon_pos,
                     default_mga=None):
    """Get MGAEntropy scores from pre-computed scores in array.

    Parameters
    ----------
    mga_vec : np.array
        numpy vector containing MGA Entropy conservation scores for residues
    codon_pos: list of int
        position of codon in protein sequence
    default_mga: float or None, default=None
        value to use if MGA entropy score not available for a given mutation.
        Drop mutations if no default specified.

    Returns
    -------
    mga_ent_scores : np.array
        score results for MGA entropy conservation
    """
    # keep only positions in range of MGAEntropy scores
    len_mga = len(mga_vec)
    good_codon_pos = [p for p in codon_pos if p < len_mga]

    # get MGAEntropy scores
    if good_codon_pos:
        mga_ent_scores = mga_vec[good_codon_pos]
    else:
        mga_ent_scores = None

    return mga_ent_scores


def read_neighbor_graph_pickle(gname, graph_dir):
    """Read in neighbor graph for given gene.

    Parameters
    ----------
    gname : str
        name of gene
    graph_dir : str
        directory containing gene graphs

    Returns
    -------
    gene_graph : dict or None
        neighbor graph as dict for gene. Returns None if not found.
    """
    graph_path = os.path.join(graph_dir, gname+".pickle")
    if os.path.exists(graph_path):
        with open(graph_path) as handle:
            gene_graph = pickle.load(handle)
        return gene_graph
    else:
        return None


def compute_ng_stat(gene_graph, pos_ct, alpha=.5):
    """Compute the clustering score for the gene on its neighbor graph.

    Parameters
    ----------
    gene_graph : dict
        Graph of spatially near codons. keys = nodes, edges = key -> value.
    pos_ct : dict
        missense mutation count for each codon
    alpha : float
        smoothing factor

    Returns
    -------
    graph_score : float
        score measuring the clustering of missense mutations in the graph
    coverage : int
        number of nodes that received non-zero weight
    """
    # skip if there are no missense mutations
    if not len(pos_ct):
        return 1.0, 0

    max_pos = max(gene_graph)
    codon_vals = np.zeros(max_pos+1)

    # smooth out mutation counts
    for pos in pos_ct:
        mut_count = pos_ct[pos]

        # update neighbor values
        neighbors = list(gene_graph[pos])
        num_neighbors = len(neighbors)
        codon_vals[neighbors] += alpha*mut_count

        # update self-value
        codon_vals[pos] += (1-alpha)*mut_count

    # compute the normalized entropy
    #total_cts = float(np.count_nonzero(codon_vals))
    #graph_score = mymath.normalized_mutation_entropy(codon_vals, total_cts=total_cts)

    # compute regular entropy
    p = codon_vals / np.sum(codon_vals)
    graph_score = mymath.shannon_entropy(p)

    # get coverage
    coverage = np.count_nonzero(p)

    return graph_score, coverage
