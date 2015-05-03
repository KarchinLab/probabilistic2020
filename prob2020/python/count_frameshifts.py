import prob2020.python.utils as utils
import prob2020.python.indel as indel
import pandas as pd


def count_frameshifts(mut_df,
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
        mut_df['End_Position'] -= 1

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

