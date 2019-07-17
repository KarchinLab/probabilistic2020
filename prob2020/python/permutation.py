import numpy as np
import csv
import prob2020.python.utils as utils
from ..cython import cutils
import prob2020.python.mutation_context as mc
import prob2020.python.scores as scores


def deleterious_permutation(obs_del,
                            context_counts,
                            context_to_mut,
                            seq_context,
                            gene_seq,
                            num_permutations=10000,
                            stop_criteria=100,
                            pseudo_count=0,
                            max_batch=25000):
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

    # calculate the # of batches for simulations
    max_batch = min(num_permutations, max_batch)
    num_batches = num_permutations // max_batch
    remainder = num_permutations % max_batch
    batch_sizes = [max_batch] * num_batches
    if remainder:
        batch_sizes += [remainder]

    num_sim = 0
    null_del_ct = 0
    for j, batch_size in enumerate(batch_sizes):
        # stop iterations if reached sufficient precision
        if null_del_ct >= stop_criteria:
            #j = j - 1
            break

        # get random positions determined by sequence context
        tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                                batch_size)
        tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

        # determine result of random positions
        for i, row in enumerate(tmp_mut_pos):
            # get info about mutations
            tmp_mut_info = mc.get_aa_mut_info(row,
                                              somatic_base,
                                              gene_seq)

            # calc deleterious mutation info
            tmp_del_count = cutils.calc_deleterious_info(tmp_mut_info['Reference AA'],
                                                         tmp_mut_info['Somatic AA'],
                                                         tmp_mut_info['Codon Pos'])

            # update empricial null distribution
            if tmp_del_count >= obs_del: null_del_ct += 1

            # stop if reach sufficient precision on p-value
            if null_del_ct >= stop_criteria:
                break
        # update number of simulations
        num_sim += i + 1

    #num_sim = j*max_batch + i+1
    del_pval = float(null_del_ct) / (num_sim)

    return del_pval


def position_permutation(obs_stat,
                         context_counts,
                         context_to_mut,
                         seq_context,
                         gene_seq,
                         gene_vest=None,
                         num_permutations=10000,
                         stop_criteria=100,
                         pseudo_count=0,
                         max_batch=25000):
    """Performs null-permutations for position-based mutation statistics
    in a single gene.

    Parameters
    ----------
    obs_stat : tuple, (recur ct, entropy, delta entropy, mean vest)
        tuple containing the observed statistics
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
    stop_criteria : int
        stop after stop_criteria iterations are more significant
        then the observed statistic.
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
    # get contexts and somatic base
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # calculate the # of batches for simulations
    max_batch = min(num_permutations, max_batch)
    num_batches = num_permutations // max_batch
    remainder = num_permutations % max_batch
    batch_sizes = [max_batch] * num_batches
    if remainder:
        batch_sizes += [remainder]

    obs_recur, obs_ent, obs_delta_ent, obs_vest = obs_stat
    num_sim = 0 # number of simulations
    null_num_recur_ct, null_entropy_ct, null_delta_entropy_ct, null_vest_ct = 0, 0, 0, 0
    for j, batch_size in enumerate(batch_sizes):
        # stop iterations if reached sufficient precision
        if null_vest_ct >= stop_criteria and null_entropy_ct >= stop_criteria:
            break

        # get random positions determined by sequence context
        tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                                batch_size)
        tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

        # calculate position-based statistics as a result of random positions
        for i, row in enumerate(tmp_mut_pos):
            # get info about mutations
            tmp_mut_info = mc.get_aa_mut_info(row,
                                              somatic_base,
                                              gene_seq)

            # calculate position info
            tmp_recur_ct, tmp_entropy, tmp_delta_entropy, _ = cutils.calc_pos_info(tmp_mut_info['Codon Pos'],
                                                                                tmp_mut_info['Reference AA'],
                                                                                tmp_mut_info['Somatic AA'],
                                                                                pseudo_count=pseudo_count,
                                                                                is_obs=0)
            # get vest scores
            if gene_vest:
                tmp_vest = scores.compute_vest_stat(gene_vest,
                                                    tmp_mut_info['Reference AA'],
                                                    tmp_mut_info['Somatic AA'],
                                                    tmp_mut_info['Codon Pos'])
            else:
                tmp_vest = 0.0

            # update empirical null distribution counts
            if tmp_entropy-utils.epsilon <= obs_ent: null_entropy_ct += 1
            if tmp_vest+utils.epsilon >= obs_vest: null_vest_ct += 1

            # stop iterations if reached sufficient precision
            if null_vest_ct >= stop_criteria and null_entropy_ct >= stop_criteria:
                break
        # update the number of simulations
        num_sim += i+1

    # calculate p-value from empirical null-distribution
    ent_pval = float(null_entropy_ct) / (num_sim)
    vest_pval = float(null_vest_ct) / (num_sim)

    return ent_pval, vest_pval


def hotmaps_permutation(obs_stat,
                        context_counts,
                        context_to_mut,
                        seq_context,
                        gene_seq,
                        window,
                        num_permutations=10000,
                        stop_criteria=100,
                        max_batch=25000,
                        null_save_path=None):
    """Performs null-permutations for position-based mutation statistics
    in a single gene.

    Parameters
    ----------
    obs_stat : dict
        dictionary mapping codons to the sum of mutations in a window
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
    window : int
        Number of codons to the left/right of a mutated position to consider
        in the window
    num_permutations : int, default: 10000
        number of permutations to create for null
    stop_criteria : int
        stop after stop_criteria iterations are more significant
        then the observed statistic.
    max_batch : int
        maximum number of whole gene simulations to do at once.
        For large number of simulations holding a matrix of M x N,
        where M is the number of mutations and N is the number of simulations,
        can get quite large.
    null_save_path : str or None
        File path to save null distribution. If None, don't save it.

    Returns
    -------
    pvals : dict
        Maps mutated codon position to the calculated p-value
    """
    # get contexts and somatic base
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # calculate the # of batches for simulations
    max_batch = min(num_permutations, max_batch)
    num_batches = num_permutations // max_batch
    remainder = num_permutations % max_batch
    batch_sizes = [max_batch] * num_batches
    if remainder:
        batch_sizes += [remainder]

    # figure out which position has highest value
    max_key = {w: max(obs_stat[w], key=(lambda key: obs_stat[w][key]))
               for w in window}

    # setup null dist counts
    null_cts = {w: {k: 0 for k in obs_stat[w]}
                for w in window }

    # empirical null distribution (saved if file path provided)
    empirical_null = {w: {} for w in window}

    num_sim = 0 # number of simulations
    for j, batch_size in enumerate(batch_sizes):
        # stop iterations if reached sufficient precision
        stop_flag = [(null_cts[w][max_key[w]]>=stop_criteria)
                      for w in window]
        if all(stop_flag):
            break
        #if null_cts[max_key] >= stop_criteria:
            #break

        # get random positions determined by sequence context
        tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                                batch_size)
        tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

        # calculate position-based statistics as a result of random positions
        for i, row in enumerate(tmp_mut_pos):
            # get info about mutations
            tmp_mut_info = mc.get_aa_mut_info(row,
                                              somatic_base,
                                              gene_seq)

            # calculate position info
            tmp_pos, tmp_sim = utils.calc_windowed_sum(tmp_mut_info['Codon Pos'],
                                                 tmp_mut_info['Reference AA'],
                                                 tmp_mut_info['Somatic AA'],
                                                 window)

            # update the counts when the empirical null passes the observed
            for tmp_w in tmp_sim:
                for tmp_key in tmp_sim[tmp_w]:
                    # get mutation count for simulation
                    val = tmp_sim[tmp_w][tmp_key]

                    # add to empirical null distribution
                    empirical_null[tmp_w].setdefault(val, 0)
                    empirical_null[tmp_w][val] += 1

                    # update counts used for p-value
                    for key in null_cts[tmp_w]:
                        if val >= obs_stat[tmp_w][key]:
                            null_cts[tmp_w][key] += 1

            # update the number of simulations
            num_sim += len(tmp_pos)

            # stop iterations if reached sufficient precision
            stop_flag = [(null_cts[w][max_key[w]]>=stop_criteria)
                         for w in window]
            if all(stop_flag):
                break

    # calculate p-value from empirical null-distribution
    pvals = {w: {k: float(null_cts[w][k]) / (num_sim) for k in obs_stat[w]}
             for w in window}

    # save empirical distribution
    if null_save_path:
        for w in window:
            # create null distribution
            output = [['mutation_count', 'p-value']]
            sorted_cts = sorted(empirical_null[w].keys())
            tmp_sum = 0
            for i in range(len(sorted_cts)):
                tmp_sum += empirical_null[w][sorted_cts[-(i+1)]]
                tmp_pval = tmp_sum / float(num_sim)
                output.append([sorted_cts[-(i+1)], tmp_pval])
            # save output
            with open(null_save_path.format(w), 'w') as handle:
                mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
                mywriter.writerows(output)

    return pvals


def protein_permutation(graph_score,
                        num_codons_obs,
                        context_counts,
                        context_to_mut,
                        seq_context,
                        gene_seq,
                        gene_graph,
                        num_permutations=10000,
                        stop_criteria=100,
                        pseudo_count=0):
    """Performs null-simulations for position-based mutation statistics
    in a single gene.

    Parameters
    ----------
    graph_score : float
        clustering score for observed data
    num_codons_obs : int
        number of codons with missense mutation in observed data
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
    stop_criteria : int
        stop after stop_criteria iterations are more significant
        then the observed statistic.

    Returns
    -------
    protein_pval : float
        p-value for clustering in neighbor graph constructure from protein
        structures
    """
    # get contexts and somatic base
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

    # calculate position-based statistics as a result of random positions
    null_graph_entropy_ct = 0
    coverage_list = []
    num_mut_list = []
    graph_entropy_list = []
    for i, row in enumerate(tmp_mut_pos):
        # calculate the expected value of the relative increase in coverage
        if i == stop_criteria-1:
            rel_inc = [coverage_list[k] / float(num_mut_list[k])
                       for k in range(stop_criteria-1)
                       if coverage_list[k]]
            exp_rel_inc = np.mean(rel_inc)

            # calculate observed statistic
            if num_codons_obs:
                obs_stat = graph_score / np.log2(exp_rel_inc*num_codons_obs)
            else:
                obs_stat = 1.0

            # calculate statistics for simulated data
            sim_stat_list = [ent / np.log2(exp_rel_inc*num_mut_list[l])
                             for l, ent in enumerate(graph_entropy_list)]
            null_graph_entropy_ct = len([s for s in sim_stat_list
                                         if s-utils.epsilon <= obs_stat])

        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # calculate position info
        tmp_tuple = cutils.calc_pos_info(tmp_mut_info['Codon Pos'],
                                         tmp_mut_info['Reference AA'],
                                         tmp_mut_info['Somatic AA'],
                                         pseudo_count=pseudo_count,
                                         is_obs=0)
        _, _, _, tmp_pos_ct = tmp_tuple

        # record num of mut codons
        if i < stop_criteria-1:
            tmp_num_mut_codons = len(tmp_pos_ct)
            num_mut_list.append(tmp_num_mut_codons)

        # get entropy on graph-smoothed probability distribution
        tmp_graph_entropy, tmp_coverage = scores.compute_ng_stat(gene_graph, tmp_pos_ct)

        # record the "coverage" in the graph
        if i < stop_criteria-1:
            coverage_list.append(tmp_coverage)
            graph_entropy_list.append(tmp_graph_entropy)

        # update empirical null distribution counts
        if i >= stop_criteria:
            #if tmp_graph_entropy-utils.epsilon <= graph_score:
            if tmp_num_mut_codons:
                sim_stat = tmp_graph_entropy / np.log2(exp_rel_inc*tmp_num_mut_codons)
            else:
                sim_stat = 1.0

            # add count
            if sim_stat-utils.epsilon <= obs_stat:
                null_graph_entropy_ct += 1

        # stop iterations if reached sufficient precision
        if null_graph_entropy_ct >= stop_criteria:
            break

    # calculate p-value from empirical null-distribution
    protein_pval = float(null_graph_entropy_ct) / (i+1)

    return protein_pval, obs_stat


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
    tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

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
    tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

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


def summary_permutation(context_counts,
                        context_to_mut,
                        seq_context,
                        gene_seq,
                        score_dir,
                        num_permutations=10000,
                        min_frac=0.0,
                        min_recur=2,
                        drop_silent=False):
    """Performs null-permutations and summarizes the results as features over
    the gene.

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
    drop_silent : bool, default=False
        Flage on whether to drop all silent mutations. Some data sources
        do not report silent mutations, and the simulations should match this.

    Returns
    -------
    summary_info_list : list of lists
        list of non-silent and silent mutation counts under the null along
        with information on recurrent missense counts and missense positional
        entropy.
    """
    mycontexts = context_counts.index.tolist()
    somatic_base = [base
                    for one_context in mycontexts
                    for base in context_to_mut[one_context]]

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

    # determine result of random positions
    gene_name = gene_seq.bed.gene_name
    gene_len = gene_seq.bed.cds_len
    summary_info_list = []
    for i, row in enumerate(tmp_mut_pos):
        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # Get all metrics summarizing each gene
        tmp_summary = cutils.calc_summary_info(tmp_mut_info['Reference AA'],
                                               tmp_mut_info['Somatic AA'],
                                               tmp_mut_info['Codon Pos'],
                                               gene_name,
                                               score_dir,
                                               min_frac=min_frac,
                                               min_recur=min_recur)

        # drop silent if needed
        if drop_silent:
            # silent mutation count is index 1
            tmp_summary[1] = 0

        # limit the precision of floats
        #pos_ent = tmp_summary[-1]
        #tmp_summary[-1] = '{0:.5f}'.format(pos_ent)

        summary_info_list.append([gene_name, i+1, gene_len]+tmp_summary)
    return summary_info_list


def maf_permutation(context_counts,
                    context_to_mut,
                    seq_context,
                    gene_seq,
                    num_permutations=10000,
                    drop_silent=False):
    """Performs null-permutations across all genes and records the results in
    a format like a MAF file. This could be useful for examining the null
    permutations because the alternative approaches always summarize the results.
    With the simulated null-permutations, novel metrics can be applied to create
    an empirical null-distribution.

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
    drop_silent : bool, default=False
        Flage on whether to drop all silent mutations. Some data sources
        do not report silent mutations, and the simulations should match this.

    Returns
    -------
    maf_list : list of tuples
        list of null mutations with mutation info in a MAF like format
    """
    mycontexts = context_counts.index.tolist()
    somatic_base, base_context = zip(*[(base, one_context)
                                       for one_context in mycontexts
                                       for base in context_to_mut[one_context]])

    # get random positions determined by sequence context
    tmp_contxt_pos = seq_context.random_pos(context_counts.iteritems(),
                                            num_permutations)
    tmp_mut_pos = np.hstack([pos_array for base, pos_array in tmp_contxt_pos])

    # info about gene
    gene_name = gene_seq.bed.gene_name
    strand = gene_seq.bed.strand
    chrom = gene_seq.bed.chrom
    gene_seq.bed.init_genome_coordinates()  # map seq pos to genome

    # determine result of random positions
    maf_list = []
    for row in tmp_mut_pos:
        # get genome coordinate
        pos2genome = np.vectorize(lambda x: gene_seq.bed.seqpos2genome[x]+1)
        genome_coord = pos2genome(row)

        # get info about mutations
        tmp_mut_info = mc.get_aa_mut_info(row,
                                          somatic_base,
                                          gene_seq)

        # get string describing variant
        var_class = cutils.get_variant_classification(tmp_mut_info['Reference AA'],
                                                      tmp_mut_info['Somatic AA'],
                                                      tmp_mut_info['Codon Pos'])

        # prepare output
        for k, mysomatic_base in enumerate(somatic_base):
            # format DNA change
            ref_nuc = tmp_mut_info['Reference Nuc'][k]
            nuc_pos = row[k]
            dna_change = 'c.{0}{1}>{2}'.format(ref_nuc, nuc_pos, mysomatic_base)

            # format protein change
            ref_aa = tmp_mut_info['Reference AA'][k]
            somatic_aa = tmp_mut_info['Somatic AA'][k]
            codon_pos = tmp_mut_info['Codon Pos'][k]
            protein_change = 'p.{0}{1}{2}'.format(ref_aa, codon_pos, somatic_aa)

            # reverse complement if on negative strand
            if strand == '-':
                ref_nuc = utils.rev_comp(ref_nuc)
                mysomatic_base = utils.rev_comp(mysomatic_base)

            # append results
            if drop_silent and var_class[k].decode() == 'Silent': continue
            maf_line = [gene_name, strand, chrom, genome_coord[k], genome_coord[k],
                        ref_nuc, mysomatic_base, base_context[k], dna_change,
                        protein_change, var_class[k].decode()]
            maf_list.append(maf_line)

    return maf_list
