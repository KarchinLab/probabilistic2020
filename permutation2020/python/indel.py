"""This module handles indels and frameshift mutations.

Indels and frameshifts are detected from the allele columns of the mutation
input.
"""
import numpy as np


def reassign_indels(indel_lens, indel_counts, gene_length):
    """Reassign indel counts to different genes based on each gene's length.

    Parameters
    ----------
    indel_lens : np.array or iterable
        length of indels to be considered
    indel_counts : np.array or iterable
        number of indels corresponding to indel_lens
    gene_length : pd.Series
        gene lengths with each gene name as an index for the pandas series

    Returns
    -------
    gene_indels : dict
        indel counts reassigned to different genes
    """
    gene_indels = {}
    prng = np.random.RandomState()
    gene_prob = gene_length / gene_length.sum()
    for i in range(len(indel_lens)):
        # randomly reassign indels based only on gene length
        gene_cts = prng.multinomial(indel_counts[i], gene_prob, size=1)
        gene_ix = np.nonzero(gene_cts)

        # record counts for indels in gene_indels dictionary
        for g in gene_ix:
            gene_name = gene_length.index[g]
            gene_indels.setdefault(gene_name, {})
            gene_indels[gene_name][indel_lens[i]] = gene_cts[g]
    return gene_indels


def keep_indels(mut_df,
                indel_len_col=True):
    """Filters out all mutations that are not indels.

    Requires that one of the alleles have '-' indicating either an insertion
    or deletion depending if found in reference allele or somatic allele
    columns, respectively.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input file as a dataframe in standard format
    indel_len_col : bool
        whether or not to add a column indicating the length of the indel

    Returns
    -------
    mut_df : pd.DataFrame
        mutations with only frameshift mutations kept
    """
    if indel_len_col:
        # calculate length
        mut_df['indel len'] = mut_df['End_Position'] - mut_df['Start_Position']

    # keep only frameshifts
    mut_df = mut_df[is_indel(mut_df)]
    return mut_df


def keep_frameshifts(mut_df,
                     indel_len_col=True):
    """Filters out all mutations that are not frameshift indels.

    Requires that one of the alleles have '-' indicating either an insertion
    or deletion depending if found in reference allele or somatic allele
    columns, respectively.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input file as a dataframe in standard format
    indel_len_col : bool
        whether or not to add a column indicating the length of the frameshift

    Returns
    -------
    mut_df : pd.DataFrame
        mutations with only frameshift mutations kept
    """
    if indel_len_col:
        # calculate length
        mut_df['indel len'] = mut_df['End_Position'] - mut_df['Start_Position']

    # keep only frameshifts
    mut_df = mut_df[is_frameshift(mut_df)]
    return mut_df


def is_frameshift(mut_df):
    """Simply returns a series indicating whether each corresponding mutation
    is a frameshift.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input file as a dataframe in standard format

    Returns
    -------
    is_fs : pd.Series
        pandas series indicating if mutaitons are frameshifts
    """
    # calculate length, 0-based coordinates
    indel_len = mut_df['End_Position'] - mut_df['Start_Position']

    # only non multiples of 3 are frameshifts
    is_fs = (indel_len%3)>0

    # make sure no single base substitutions are counted
    is_indel = (mut_df['Reference_Allele']=='-') | (mut_df['Tumor_Allele']=='-')
    is_fs[~is_indel] = False
    return is_fs


def is_indel(mut_df):
    """Simply returns a series indicating whether each corresponding mutation
    is an indel.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input file as a dataframe in standard format

    Returns
    -------
    is_indel : pd.Series
        pandas series indicating if mutaitons are indels
    """
    # calculate length, 0-based coordinates
    indel_len = mut_df['End_Position'] - mut_df['Start_Position']

    # make sure no single base substitutions are counted
    is_indel = (mut_df['Reference_Allele']=='-') | (mut_df['Tumor_Allele']=='-')

    # make sure indel has a length
    is_indel[indel_len<1] = False
    return is_indel


def get_frameshift_lengths(num_bins):
    """Simple function that returns the lengths for each frameshift category
    if `num_bins` number of frameshift categories are requested.
    """
    fs_len = []
    i = 1
    tmp_bins = 0
    while(tmp_bins<num_bins):
        if i%3:
            fs_len.append(i)
            tmp_bins += 1
        i += 1
    return fs_len
