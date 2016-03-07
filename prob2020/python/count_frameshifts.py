import prob2020.python.utils as utils
import prob2020.python.indel as indel
import pandas as pd


def count_frameshift_total(mut_df,
                           bed_path,
                           use_unmapped=False,
                           to_zero_based=False):
    """Count frameshifts for each gene.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input
    bed_path : str
        path to BED file containing reference tx for genes
    use_unmapped : Bool
        flag indicating whether to include frameshifts not mapping
        to reference tx
    to_zero_based : Bool
        whether to convert end-coordinate to zero based for analysis

    Returns
    -------
    fs_cts_df : pd.DataFrame
        contains both total frameshift counts and frameshift counts
        not mappable to the reference transcript.
    """
    if to_zero_based:
        mut_df['Start_Position'] = mut_df['Start_Position'] - 1

    fs_cts = {}  # frameshift count information for each gene
    fs_df = indel.keep_frameshifts(mut_df)

    for bed in utils.bed_generator(bed_path):
        gene_df = fs_df[fs_df['Gene']==bed.gene_name]

        # find it frameshift actually is on gene annotation
        fs_pos = []
        for ix, row in gene_df.iterrows():
            indel_pos = [row['Start_Position'], row['End_Position']]
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], indel_pos)
            fs_pos.append(coding_pos)

        # mark frameshifts that could not be mapped to reference tx
        gene_df['unmapped'] = [(1 if x is None else 0) for x in fs_pos]
        total_fs = len(gene_df)
        unmapped_fs = len(gene_df[gene_df['unmapped']==1])

        # filter out frameshifts that did not match reference tx
        if not use_unmapped:
            gene_df = gene_df[gene_df['unmapped']==0]
            total_fs -= unmapped_fs

        info = [total_fs, unmapped_fs,]
        fs_cts[bed.gene_name] = info

    # prepare counts into a dataframe
    fs_cts_df = pd.DataFrame.from_dict(fs_cts, orient='index')
    cols = ['total', 'unmapped',]
    fs_cts_df.columns = cols

    return fs_cts_df


def count_frameshift_bins(mut_df,
                          bed_path,
                          num_bins,
                          num_samples=None,
                          use_unmapped=False,
                          to_zero_based=False):
    """Count frameshifts for each gene.

    Parameters
    ----------
    mut_df : pd.DataFrame
        mutation input
    bed_path : str
        path to BED file containing reference tx for genes
    num_bins : int
        number of bins to stratify frameshift lengths
    num_samples : int
        number of samples where mutations were gathered
    use_unmapped : Bool
        flag indicating whether to include frameshifts not mapping
        to reference tx
    to_zero_based : Bool
        whether to convert end-coordinate to zero based for analysis

    Returns
    -------
    fs_cts_df : pd.DataFrame
        contains both total frameshift counts and frameshift counts
        not mappable to the reference transcript.
    """
    # convert to zero-based
    if to_zero_based:
        mut_df['Start_Position'] = mut_df['Start_Position'] - 1

    # number of samples
    if num_samples is None:
        num_samples = mut_df['Tumor_Sample'].nunique()

    fs_cts = {}  # frameshift count information for each gene
    fs_df = indel.keep_frameshifts(mut_df)
    fs_lens = indel.get_frameshift_lengths(num_bins)

    for bed in utils.bed_generator(bed_path):
        gene_df = fs_df[fs_df['Gene']==bed.gene_name]

        # find it frameshift actually is on gene annotation
        fs_pos = []
        for ix, row in gene_df.iterrows():
            indel_pos = [row['Start_Position'], row['End_Position']]
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], indel_pos)
            fs_pos.append(coding_pos)

        # mark frameshifts that could not be mapped to reference tx
        gene_df['unmapped'] = [(1 if x is None else 0) for x in fs_pos]
        total_fs = len(gene_df)
        unmapped_fs = len(gene_df[gene_df['unmapped']==1])

        # filter out frameshifts that did not match reference tx
        if not use_unmapped:
            gene_df = gene_df[gene_df['unmapped']==0]

        # count all frameshifts
        all_len_counts = gene_df['indel len'].value_counts()
        all_lens = list(set(all_len_counts.index) | set(fs_lens))
        all_len_counts = all_len_counts.reindex(all_lens).fillna(0)
        # count only in specified lengths
        fs_len_cts = all_len_counts.reindex(fs_lens)
        fs_len_cts[max(fs_len_cts.index)] = all_len_counts[all_len_counts.index>=fs_lens[-1]].sum()
        fs_len_cts = fs_len_cts.fillna(0).astype(int)

        # get length of gene
        gene_len = bed.cds_len
        gene_bases_at_risk = gene_len * num_samples

        info = [total_fs, unmapped_fs, gene_len, gene_bases_at_risk]
        fs_cts[bed.gene_name] = fs_len_cts.tolist() + info

    fs_cts_df = pd.DataFrame.from_dict(fs_cts, orient='index')
    cols = list(map(str, fs_lens)) + ['total', 'unmapped', 'gene length', 'bases at risk']
    fs_cts_df.columns = cols
    return fs_cts_df
