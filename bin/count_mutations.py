#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import prob2020.python.utils as utils
import prob2020.python.mutation_context as mc
from prob2020.python.gene_sequence import GeneSequence
import prob2020.cython.cutils as cutils
import pandas as pd
import pysam
import argparse


def keep_inframe_indels(mut_df):
    # calculate length
    mut_df['indel len'] = mut_df['End_Position'] - mut_df['Start_Position']

    # filter out non-indels
    indel_ixs = (mut_df['Reference_Allele']=='-') | (mut_df['Tumor_Allele']=='-')

    # filter out non-frameshift indels
    indel_df = mut_df[(indel_ixs) & ((mut_df['indel len']%3)==0)]

    # drop all indels after extracting inframe indels
    mut_df = mut_df[~indel_ixs]

    return mut_df, indel_df


def count_mutations(mut_df,
                    bed_path,
                    num_samples,
                    opts,
                    use_unmapped=False):
    mut_cts = {}  # count information for each gene
    mut_df, indel_df = keep_inframe_indels(mut_df)
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])
    all_contexts = mc.get_all_context_names(opts['context'])

    # pre-get the dataframe split into chunks by gene name
    gene_grp = mut_df.groupby('Gene')
    gene_mut_data = {mygene: mut_df.ix[gene_grp.groups[mygene], :]
                     for mygene in gene_grp.groups}

    for bed in utils.bed_generator(bed_path):
        # get mutation data for gene
        if bed.gene_name in gene_mut_data:
            # there are mutations for this gene case
            gene_mut_df = gene_mut_data[bed.gene_name]
        else:
            # no data case
            cols = ['Gene', 'Chromosome', 'Start_Position', 'Reference_Allele',
                    'Tumor_Allele', 'Variant_Classification', 'Protein_Change']
            gene_mut_df = pd.DataFrame(columns=cols)
        gene_indel_df = indel_df[indel_df['Gene']==bed.gene_name]
        gs.set_gene(bed)

        # get info on SNV counts
        mut_info = mc.compute_mutation_context(bed, gs, mut_df, opts)
        context_cts, context_to_mut, mut_base, gs, sc  = mut_info
        seq_context_bases = [(len(sc.context2pos[ctxt]) if sc.is_valid_context(ctxt) else 0)
                             for ctxt in all_contexts]
        seq_context_bases_at_risk = [c*num_samples for c in seq_context_bases]
        all_context_counts = []
        for context in all_contexts:
            tmp_mut = mut_base[mut_base['Context']==context]
            if len(tmp_mut) > 0:
                aa_info = mc.get_aa_mut_info(tmp_mut['Coding Position'].tolist(),
                                             tmp_mut['Tumor_Allele'].tolist(),
                                             gs)
                tmp_non_silent = cutils.calc_non_silent_info(aa_info['Reference AA'],
                                                             aa_info['Somatic AA'])
                num_non_silent = tmp_non_silent[0]
            else:
                num_non_silent = 0
            all_context_counts.append(num_non_silent)

        # find position of indels and SNVs on gene
        indel_pos_list = []
        mut_pos_list = []
        for ix, row in gene_indel_df.iterrows():
            # add indels
            indel_pos = [row['Start_Position'], row['End_Position']]
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], indel_pos)
            indel_pos_list.append(coding_pos)

        # mark SNVs and indels that could not be mapped to reference tx
        gene_indel_df['unmapped'] = [(1 if x is None else 0)
                                     for x in indel_pos_list]
        unmapped_indel = len(gene_indel_df[gene_indel_df['unmapped']==1])

        # filter out mutations/indels that did not match reference tx
        if not use_unmapped:
            total_indel = len(gene_indel_df[gene_indel_df['unmapped']==0])
        else:
            total_indel = len(gene_indel_df)

        # get length of gene
        gene_len = bed.cds_len
        indel_bases_at_risk = gene_len * num_samples

        indel_info = [total_indel, unmapped_indel,
                      gene_len, indel_bases_at_risk]
        mut_cts[bed.gene_name] = all_context_counts + seq_context_bases_at_risk + indel_info
    gene_fa.close()  # close fasta

    # format counts into a dataframe
    mut_cts_df = pd.DataFrame.from_dict(mut_cts, orient='index')
    mut_bases_at_risk_cols = [c + ' bases at risk' for c in all_contexts]
    indel_cols = ['indel', 'indel unmapped',
                  'gene length', 'indel bases at risk']
    cols = all_contexts + mut_bases_at_risk_cols + indel_cols
    mut_cts_df.columns = cols
    return mut_cts_df


def parse_arguments():
    info = 'Counts the number of SNVs and inframe indels in each gene'
    parser = argparse.ArgumentParser(description=info)
    help_str = 'Mutation file containing SNVs and inframe indels'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Gene annotation in BED format with a single reference transcript'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = 'gene FASTA file from extract_gene_seq.py script'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Number of DNA bases to use as context. 0 indicates no context. '
                '1 indicates only use the mutated base.  1.5 indicates using '
                'the base context used in CHASM '
                '(http://wiki.chasmsoftware.org/index.php/CHASM_Overview). '
                '2 indicates using the mutated base and the upstream base. '
                '3 indicates using the mutated base and both the upstream '
                'and downstream bases. (Default: 1.5)')
    parser.add_argument('-c', '--context',
                        type=float, default=1.5,
                        help=help_str)
    help_str = 'Number of sequenced samples'
    parser.add_argument('-n', '--sample-number',
                        type=int, required=True,
                        help=help_str)
    help_str = 'Use mutations that could not be placed onto the reference transcript'
    parser.add_argument('-u', '--use-unmapped',
                        action='store_true', default=False,
                        help=help_str)
    help_str = 'Output of mutation counts'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    df = pd.read_csv(opts['mutations'], sep='\t')
    df['Start_Position'] = df['Start_Position'] - 1  # convert to 0-based coord

    # count mutations
    mut_cts = count_mutations(df, opts['bed'],
                              opts['sample_number'],
                              opts, opts['use_unmapped'])

    mut_cts = mut_cts.drop(['indel', 'indel unmapped',
                            'gene length', 'indel bases at risk'],
                            axis=1)

    # save results
    if opts['output']:
        mut_cts.to_csv(opts['output'], sep='\t')

    return mut_cts


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
