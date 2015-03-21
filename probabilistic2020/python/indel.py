"""This module handles indels and frameshift mutations.

Indels and frameshifts are detected from the allele columns of the mutation
input.
"""
import numpy as np
import pandas as pd

def simulate_indel_counts(indel_df, bed_dict, num_permutations=1):
    # count indels
    bed_genes = [mybed
                 for chrom in bed_dict
                 for mybed in bed_dict[chrom]]
    tmp = []
    for b in bed_genes:
        b.init_genome_coordinates()
        tmp.append(b)
    bed_genes = tmp
    gene_lengths = pd.Series([b.cds_len for b in bed_genes],
                              index=[b.gene_name for b in bed_genes])

    # generate random indel assignments
    gene_prob = gene_lengths.astype(float) / gene_lengths.sum()
    indel_lens = indel_df['indel len'].copy().values
    is_fs = (indel_lens % 3) > 0
    indel_ixs = np.arange(len(indel_lens))
    prng = np.random.RandomState(seed=None)

    # randomly reassign indels
    mygene_cts = prng.multinomial(len(indel_lens), gene_prob, size=num_permutations)
    inframe_cts = mygene_cts.copy()
    for row in range(mygene_cts.shape[0]):
        nonzero_ix = np.nonzero(mygene_cts[row,:])[0]

        # randomly shuffle indel lengths
        prng.shuffle(indel_ixs)
        is_fs = is_fs[indel_ixs]

        # iterate over each gene
        indel_ix = 0
        for j in range(len(nonzero_ix)):
            prev_indel_ix = indel_ix
            num_gene_indels = mygene_cts[row, nonzero_ix[j]]
            indel_ix += num_gene_indels
            inframe_cts[row, nonzero_ix[j]] = num_gene_indels - np.sum(is_fs[prev_indel_ix:indel_ix])
    mygene_cts -= inframe_cts
    return mygene_cts, inframe_cts, gene_lengths.index


def simulate_indel_maf(indel_df, bed_dict, num_permutations=1):
    # count indels
    bed_genes = [mybed
                for chrom in bed_dict
                for mybed in bed_dict[chrom]]
    tmp = []
    for b in bed_genes:
        b.init_genome_coordinates()
        tmp.append(b)
    bed_genes = tmp
    gene_lengths = pd.Series([b.cds_len for b in bed_genes],
                                index=[b.gene_name for b in bed_genes])

    # generate random indel assignments
    gene_prob = gene_lengths.astype(float) / gene_lengths.sum()
    indel_lens = indel_df['indel len'].copy().values
    indel_types = indel_df['indel type'].copy().values
    indel_ixs = np.arange(len(indel_lens))
    prng = np.random.RandomState(seed=None)

    for i in range(num_permutations):
        # randomly reassign indels
        mygene_cts = prng.multinomial(len(indel_lens), gene_prob)
        nonzero_ix = np.nonzero(mygene_cts)[0]

        # randomly shuffle indel lengths
        prng.shuffle(indel_ixs)
        indel_lens = indel_lens[indel_ixs]
        indel_types = indel_types[indel_ixs]

        # iterate over each gene
        indel_ix = 0
        for j in range(len(nonzero_ix)):
            prev_indel_ix = indel_ix
            num_gene_indels = mygene_cts[nonzero_ix[j]]
            indel_ix += num_gene_indels
            maf_lines = counts2maf(num_gene_indels,
                                   indel_lens[prev_indel_ix:indel_ix],
                                   indel_types[prev_indel_ix:indel_ix],
                                   bed_genes[nonzero_ix[j]])
            yield maf_lines


def counts2maf(num_indels, myindel_lens, myindel_types, gene_bed, seed=None):
    maf_list = []
    prng = np.random.RandomState(seed=seed)
    pos = prng.randint(low=0, high=gene_bed.cds_len, size=num_indels)
    genome_pos = [gene_bed.seqpos2genome[p] for p in pos]
    is_frame_shift = myindel_lens%3
    for i, gpos in enumerate(genome_pos):
        if myindel_types[i] == 'INS':
            var_class = 'Frame_Shift_Ins' if is_frame_shift[i] else 'In_Frame_Ins'
            dna_change = 'c.{0}_{1}ins'.format(pos[i], pos[i])
            prot_change = 'p.?'
            tmp = [gene_bed.gene_name, gene_bed.strand, gene_bed.chrom,
                   gpos, gpos, '-', 'N'*myindel_lens[i], '-', dna_change,
                   prot_change, var_class]
            maf_list.append(tmp)
        else:
            var_class = 'Frame_Shift_Del' if is_frame_shift[i] else 'In_Frame_Del'
            dna_change = 'c.{0}_{1}del'.format(pos[i], pos[i]+myindel_lens[i])
            prot_change = 'p.?'
            tmp = [gene_bed.gene_name, gene_bed.strand, gene_bed.chrom,
                   gpos+1, gpos+myindel_lens[i], 'N'*myindel_lens[i], '-', '-', dna_change,
                   prot_change, var_class]
            maf_list.append(tmp)

    return maf_list


def keep_indels(mut_df,
                indel_len_col=True,
                indel_type_col=True):
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

    if indel_type_col:
        is_ins = mut_df['Reference_Allele']=='-'
        is_del = mut_df['Tumor_Allele']=='-'
        mut_df['indel type'] = ''
        mut_df['indel type'][is_ins] = 'INS'
        mut_df['indel type'][is_del] = 'DEL'

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
    # is_indel[indel_len<1] = False
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
