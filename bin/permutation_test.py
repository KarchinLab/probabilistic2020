#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# package imports
import permutation2020.python.utils as utils
from permutation2020.python.gene_sequence import GeneSequence
from permutation2020.python.sequence_context import SequenceContext
import permutation2020.cython.cutils as cutils
import permutation2020.python.permutation as pm
import permutation2020.python.mutation_context as mc

import argparse
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool
import logging
import datetime

logger = logging.getLogger(__name__)  # module logger

def calc_deleterious_p_value(mut_info,
                             unmapped_mut_info,
                             sc,
                             gs,
                             bed,
                             num_permutations,
                             del_threshold,
                             pseudo_count):
    """Calculates the p-value for the number of inactivating SNV mutations.

    Calculates p-value based on how many permutations exceed the observed value.

    Parameters
    ----------
    mut_info : dict
        contains codon and amino acid residue information for mutations mappable
        to provided reference tx.
    unmapped_mut_info : dict
        contains codon/amino acid residue info for mutations that are NOT mappable
        to provided reference tx.
    sc : SequenceContext
        object contains the nucleotide contexts for a gene such that new random
        positions can be obtained while respecting nucleotide context.
    gs : GeneSequence
        contains gene sequence
    bed : BedLine
        just used to return gene name
    num_permutations : int
        number of permutations to perform to estimate p-value. more permutations
        means more precision on the p-value.
    """
    if len(mut_info) > 0:
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        mut_info['Context'] = mut_info['Coding Position'].apply(lambda x: sc.pos2context[x])

        # group mutations by context
        cols = ['Context', 'Tumor_Allele']
        unmapped_mut_df = pd.DataFrame(unmapped_mut_info)
        tmp_df = pd.concat([mut_info[cols], unmapped_mut_df[cols]])
        context_cts = tmp_df['Context'].value_counts()
        context_to_mutations = dict((name, group['Tumor_Allele'])
                                    for name, group in tmp_df.groupby('Context'))

        # get deleterious info for actual mutations
        aa_mut_info = mc.get_aa_mut_info(mut_info['Coding Position'],
                                         mut_info['Tumor_Allele'].tolist(),
                                         gs)
        ref_aa = aa_mut_info['Reference AA'] + unmapped_mut_info['Reference AA']
        somatic_aa = aa_mut_info['Somatic AA'] + unmapped_mut_info['Somatic AA']
        codon_pos = aa_mut_info['Codon Pos'] + unmapped_mut_info['Codon Pos']
        num_del = cutils.calc_deleterious_info(ref_aa, somatic_aa, codon_pos)

        # skip permutation test if number of deleterious mutations is not at
        # least meet some user-specified threshold
        if num_del >= del_threshold:
            # perform permutations
            null_del_list = pm.deleterious_permutation(context_cts,
                                                       context_to_mutations,
                                                       sc,  # sequence context obj
                                                       gs,  # gene sequence obj
                                                       num_permutations,
                                                       pseudo_count)

            # calculate p-value
            del_num_nulls = sum([1 for d in null_del_list
                                 if d+utils.epsilon >= num_del])
            del_p_value = del_num_nulls / float(num_permutations)
        else:
            del_p_value = None
    else:
        num_del = 0
        del_p_value = None

    result = [bed.gene_name, num_del, del_p_value]
    return result


def calc_position_p_value(mut_info,
                          unmapped_mut_info,
                          sc,
                          gs,
                          bed,
                          num_permutations,
                          pseudo_count,
                          min_recurrent,
                          min_fraction):
    if len(mut_info) > 0:
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        mut_info['Context'] = mut_info['Coding Position'].apply(lambda x: sc.pos2context[x])

        # group mutations by context
        cols = ['Context', 'Tumor_Allele']
        unmapped_mut_df = pd.DataFrame(unmapped_mut_info)
        tmp_df = pd.concat([mut_info[cols], unmapped_mut_df[cols]])
        context_cts = tmp_df['Context'].value_counts()
        context_to_mutations = dict((name, group['Tumor_Allele'])
                                    for name, group in tmp_df.groupby('Context'))

        # perform permutations
        permutation_result = pm.position_permutation(context_cts,
                                                     context_to_mutations,
                                                     sc,  # sequence context obj
                                                     gs,  # gene sequence obj
                                                     num_permutations,
                                                     pseudo_count)
        num_recur_list, pos_entropy_list, delta_pos_entropy_list = permutation_result  # unpack results

        # get recurrent info for actual mutations
        aa_mut_info = mc.get_aa_mut_info(mut_info['Coding Position'],
                                         mut_info['Tumor_Allele'].tolist(),
                                         gs)
        codon_pos = aa_mut_info['Codon Pos'] + unmapped_mut_info['Codon Pos']
        ref_aa = aa_mut_info['Reference AA'] + unmapped_mut_info['Reference AA']
        somatic_aa = aa_mut_info['Somatic AA'] + unmapped_mut_info['Somatic AA']
        num_recurrent, pos_ent, delta_pos_ent = cutils.calc_pos_info(codon_pos,
                                                                     ref_aa,
                                                                     somatic_aa,
                                                                     min_frac=min_fraction,
                                                                     min_recur=min_recurrent)

        # calculate permutation p-value
        recur_num_nulls = sum([1 for null_recur in num_recur_list
                               if null_recur+utils.epsilon >= num_recurrent])
        entropy_num_nulls = sum([1 for null_ent in pos_entropy_list
                                 if null_ent-utils.epsilon <= pos_ent])
        delta_entropy_num_nulls = sum([1 for null_ent in delta_pos_entropy_list
                                       if null_ent+utils.epsilon >= delta_pos_ent])
        recur_p_value = recur_num_nulls / float(num_permutations)
        ent_p_value = entropy_num_nulls / float(num_permutations)
        delta_ent_p_value = delta_entropy_num_nulls / float(num_permutations)
    else:
        num_recurrent = 0
        pos_ent = 0
        delta_pos_ent = 0
        recur_p_value = 1.0
        ent_p_value = 1.0
        delta_ent_p_value = 1.0
    result = [bed.gene_name, num_recurrent, pos_ent, delta_pos_ent,
              recur_p_value, ent_p_value, delta_ent_p_value]
    return result


def calc_effect_p_value(mut_info,
                        unmapped_mut_info,
                        sc,
                        gs,
                        bed,
                        num_permutations,
                        pseudo_count,
                        min_recurrent,
                        min_fraction):
    if len(mut_info) > 0:
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        mut_info['Context'] = mut_info['Coding Position'].apply(lambda x: sc.pos2context[x])

        # group mutations by context
        cols = ['Context', 'Tumor_Allele']
        unmapped_mut_df = pd.DataFrame(unmapped_mut_info)
        tmp_df = pd.concat([mut_info[cols], unmapped_mut_df[cols]])
        context_cts = tmp_df['Context'].value_counts()
        context_to_mutations = dict((name, group['Tumor_Allele'])
                                    for name, group in tmp_df.groupby('Context'))

        # perform permutations
        permutation_result = pm.effect_permutation(context_cts,
                                                   context_to_mutations,
                                                   sc,  # sequence context obj
                                                   gs,  # gene sequence obj
                                                   num_permutations,
                                                   pseudo_count)
        effect_entropy_list, recur_list, inactivating_list = permutation_result  # unpack results

        # get effect info for actual mutations
        aa_mut_info = mc.get_aa_mut_info(mut_info['Coding Position'],
                                         mut_info['Tumor_Allele'].tolist(),
                                         gs)
        codon_pos = aa_mut_info['Codon Pos'] + unmapped_mut_info['Codon Pos']
        ref_aa = aa_mut_info['Reference AA'] + unmapped_mut_info['Reference AA']
        somatic_aa = aa_mut_info['Somatic AA'] + unmapped_mut_info['Somatic AA']
        effect_ent, num_recur, num_inactivating = cutils.calc_effect_info(codon_pos,
                                                                          ref_aa,
                                                                          somatic_aa,
                                                                          min_frac=min_fraction,
                                                                          min_recur=min_recurrent)

        # calculate permutation p-value
        entropy_num_nulls = sum([1 for null_ent in effect_entropy_list
                                 if null_ent-utils.epsilon <= effect_ent])
        ent_p_value = entropy_num_nulls / float(num_permutations)
    else:
        num_recur = 0
        num_inactivating = 0
        effect_ent = 0
        ent_p_value = 1.0
    result = [bed.gene_name, num_recur, num_inactivating,
              effect_ent, ent_p_value]
    return result


@utils.log_error_decorator
def singleprocess_permutation(info):
    bed_list, mut_df, opts = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    num_permutations = opts['num_permutations']
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    result = []
    for bed in bed_list:
        # prepare info for running permutation test
        gene_mut = mut_df[mut_df['Gene']==bed.gene_name]
        cols = ['Chromosome', 'Start_Position', 'Reference_Allele',
                'Tumor_Allele', 'Variant_Classification', 'Protein_Change']
        mut_info = gene_mut[cols]
        gs.set_gene(bed)
        sc = SequenceContext(gs, seed=opts['seed'])

        # count total mutations in gene
        total_mut = len(mut_info)

        # fix nucleotide letter if gene is on - strand
        if bed.strand == '-':
            mut_info['Tumor_Allele'] = mut_info['Tumor_Allele'].map(lambda x: utils.rev_comp(x))

        # get coding positions, mutations unmapped to the reference tx will have
        # NA for a coding position
        pos_list = []
        for ix, row in mut_info.iterrows():
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], row['Start_Position'])
            pos_list.append(coding_pos)
        mut_info['Coding Position'] = pos_list

        # recover mutations that could not be mapped to the reference transcript
        # for a gene before being dropped (next step)
        unmapped_mut_info = mc.recover_unmapped_mut_info(mut_info, bed, sc, opts)

        # drop mutations wich do not map to reference tx
        mut_info = mut_info.dropna(subset=['Coding Position'])  # mutations need to map to tx
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        unmapped_muts = total_mut - len(mut_info)

        # construct sequence context
        #gs.add_germline_variants(mut_info['Reference_Allele'].tolist(),
        #                         mut_info['Coding Position'].tolist())

        # calculate results of permutation test
        if opts['kind'] == 'oncogene':
            # calculate position based permutation results
            tmp_result = calc_position_p_value(mut_info, unmapped_mut_info, sc,
                                               gs, bed, num_permutations,
                                               opts['recurrent_pseudo_count'],
                                               opts['recurrent'],
                                               opts['fraction'])
            result.append(tmp_result + [total_mut, unmapped_muts])
        elif opts['kind'] == 'tsg':
            # calculate results for deleterious mutation permutation test
            tmp_result = calc_deleterious_p_value(mut_info, unmapped_mut_info,
                                                  sc, gs, bed, num_permutations,
                                                  opts['deleterious'],
                                                  opts['deleterious_pseudo_count'])
            result.append(tmp_result + [total_mut, unmapped_muts])
        else:
            # calc results for entropy-on-effect permutation test
            tmp_result = calc_effect_p_value(mut_info, unmapped_mut_info,
                                             sc, gs, bed, num_permutations,
                                             opts['recurrent_pseudo_count'],
                                             opts['recurrent'],
                                             opts['fraction'])
            result.append(tmp_result + [total_mut, unmapped_muts])

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def multiprocess_permutation(bed_dict, mut_df, opts):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys(), key=lambda x: len(bed_dict[x]), reverse=True)
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    result_list = []
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for chrom_result in process_results:
                    result_list += chrom_result
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            info = (bed_dict[chroms[i]], mut_df, opts)
            result_list += singleprocess_permutation(info)

    return result_list


def handle_tsg_results(permutation_result):
    permutation_df = pd.DataFrame(sorted(permutation_result, key=lambda x: x[2]),
                                  columns=['gene', 'num deleterious', 'deleterious p-value',
                                           'Total Mutations', 'Unmapped to Ref Tx'])
    permutation_df['deleterious p-value'] = permutation_df['deleterious p-value'].astype('float')
    tmp_df = permutation_df[permutation_df['deleterious p-value'].notnull()]

    # get benjamani hochberg adjusted p-values
    permutation_df['deleterious BH q-value'] = np.nan
    permutation_df['deleterious BH q-value'][tmp_df.index] = utils.bh_fdr(tmp_df['deleterious p-value'])

    # sort output by p-value. due to no option to specify NaN order in
    # sort, the df needs to sorted descendingly and then flipped
    permutation_df = permutation_df.sort(columns='deleterious p-value', ascending=False)
    permutation_df = permutation_df.reindex(index=permutation_df.index[::-1])

    # order result
    permutation_df = permutation_df.set_index('gene', drop=False)
    col_order  = ['gene', 'Total Mutations', 'Unmapped to Ref Tx',
                  'num deleterious', 'deleterious p-value',
                  'deleterious BH q-value']
    return permutation_df[col_order]


def handle_oncogene_results(permutation_result, non_tested_genes):
    mycols = ['gene', 'num recurrent', 'position entropy',
              'delta position entropy',
              # 'kde position entropy', 'kde bandwidth',
              'recurrent p-value', 'entropy p-value',
              'delta entropy p-value',
              # 'kde entropy p-value', 'kde bandwidth p-value',
              'Total Mutations', 'Unmapped to Ref Tx']
    permutation_df = pd.DataFrame(sorted(permutation_result, key=lambda x: x[4]),
                                  columns=mycols)

    # get benjamani hochberg adjusted p-values
    permutation_df['recurrent BH q-value'] = utils.bh_fdr(permutation_df['recurrent p-value'])
    permutation_df['entropy BH q-value'] = utils.bh_fdr(permutation_df['entropy p-value'])
    permutation_df['delta entropy BH q-value'] = utils.bh_fdr(permutation_df['delta entropy p-value'])
    #permutation_df['kde entropy BH q-value'] = utils.bh_fdr(permutation_df['kde entropy p-value'])
    #permutation_df['kde bandwidth BH q-value'] = utils.bh_fdr(permutation_df['kde bandwidth p-value'])

    # include non-tested genes in the result
    no_test_df = pd.DataFrame(index=range(len(non_tested_genes)))
    no_test_df['Performed Recurrency Test'] = 0
    no_test_df['gene'] = non_tested_genes
    permutation_df = pd.concat([permutation_df, no_test_df])
    permutation_df['Performed Recurrency Test'] = permutation_df['Performed Recurrency Test'].fillna(1).astype(int)

    # order output
    permutation_df = permutation_df.set_index('gene', drop=False)  # make sure genes are indices
    permutation_df['num recurrent'] = permutation_df['num recurrent'].fillna(-1).astype(int)  # fix dtype isssue
    col_order = ['gene', 'Total Mutations', 'Unmapped to Ref Tx',
                 'num recurrent', 'position entropy',
                 'delta position entropy',
                 #'kde position entropy', 'kde bandwidth',
                 'recurrent p-value', 'recurrent BH q-value', 'entropy p-value', 'entropy BH q-value',
                 #'kde entropy p-value', 'kde entropy BH q-value', 'kde bandwidth p-value', 'kde bandwidth BH q-value',
                 'delta entropy p-value', 'delta entropy BH q-value',
                 'Performed Recurrency Test']
    return permutation_df[col_order]


def handle_effect_results(permutation_result):
    mycols = ['gene', 'num recurrent', 'num inactivating', 'entropy-on-effect',
              'entropy-on-effect p-value',
              'Total Mutations', 'Unmapped to Ref Tx']
    permutation_df = pd.DataFrame(sorted(permutation_result, key=lambda x: x[4]),
                                  columns=mycols)

    # get benjamani hochberg adjusted p-values
    permutation_df['entropy-on-effect BH q-value'] = utils.bh_fdr(permutation_df['entropy-on-effect p-value'])

    # order output
    permutation_df = permutation_df.set_index('gene', drop=False)  # make sure genes are indices
    permutation_df['num recurrent'] = permutation_df['num recurrent'].fillna(-1).astype(int)  # fix dtype isssue
    col_order = ['gene', 'Total Mutations', 'Unmapped to Ref Tx',
                 'num recurrent', 'num inactivating', 'entropy-on-effect',
                 'entropy-on-effect p-value', 'entropy-on-effect BH q-value']
    return permutation_df[col_order]


def parse_arguments():
    # make a parser
    info = 'Performs a permutation test on the oncogene and TSG score'
    parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')

    # program arguments
    help_str = 'gene FASTA file from extract_gene_seq.py script'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'DNA mutations file'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'BED file annotation of genes'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Number of processes to use. 0 indicates using a single '
                'process without using a multiprocessing pool '
                '(more means Faster, default: 0).')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Number of permutations for null model. p-value precision '
                'increases with more permutations, however this will also '
                'increase the run time (Default: 10000).')
    parser.add_argument('-n', '--num-permutations',
                        type=int, default=10000,
                        help=help_str)
    help_str = ('Kind of permutation test to perform ("oncogene" or "tsg"). "position-based" permutation '
                'test is intended to find oncogenes using position based statistics. '
                'The "deleterious" permutation test is intended to find tumor '
                'suppressor genes. (Default: oncogene)')
    parser.add_argument('-k', '--kind',
                        type=str, default='oncogene',
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
    help_str = ('Use mutations that are not mapped to the the single reference '
                'transcript for a gene specified in the bed file indicated by '
                'the -b option.')
    parser.add_argument('-u', '--use-unmapped',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Path to the genome fasta file. Required if --use-unmapped flag '
                'is used. (Default: None)')
    parser.add_argument('-g', '--genome',
                        type=str, default='',
                        help=help_str)
    help_str = ('Minimum number of mutations at a position for it to be '
                'considered a recurrently mutated position (Default: 3).')
    parser.add_argument('-r', '--recurrent',
                        type=int, default=3,
                        help=help_str)
    help_str = ('Fraction of total mutations in a gene. This define the '
                'minimumm number of mutations for a position to be defined '
                'as recurrently mutated (Defaul: .02).')
    parser.add_argument('-f', '--fraction',
                        type=float, default=.02,
                        help=help_str)
    help_str = ('Perform tsg permutation test if gene has '
                'at least a user specified number of deleterious mutations (default: 2)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=2,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. Values greater than one indicate all '
                'genes will be tested (Default: 1.01).')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=1.01,
                        help=help_str)
    help_str = ('Deleterious mutation pseudo-count for null distribution '
                'statistics. (Default: 0)')
    parser.add_argument('-dp', '--deleterious-pseudo-count',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Recurrent missense mutation pseudo-count for null distribution '
                'statistics. (Default: 0)')
    parser.add_argument('-rp', '--recurrent-pseudo-count',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Specify the seed for the pseudo random number generator. '
                'By default, the seed is randomly chosen based. The seed will '
                'be used for the permutation test monte carlo simulations.')
    parser.add_argument('-seed', '--seed',
                        type=int, default=None,
                        help=help_str)
    help_str = 'Output of probabilistic 20/20 results'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    utils.start_logging(log_file=log_file,
                        log_level=log_level)  # start logging

    opts = vars(args)
    if opts['use_unmapped'] and not opts['genome']:
        print('You must specify a genome fasta with -g if you set the '
              '--use-unmapped flag to true.')
        sys.exit(1)

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts, mut_df=None):
    # hack to index the FASTA file
    gene_fa = pysam.Fastafile(opts['input'])
    gene_fa.close()

    # Get Mutations
    if mut_df is None:
        mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # specify genes to skip
    if opts['kind'] == 'oncogene':
        # find genes with tsg score above threshold to filter out for oncogene
        # permutation test
        non_tested_genes = utils._get_high_tsg_score(mut_df, opts['tsg_score'])
    else:
        # don't filter out genes for tsg permutation test
        non_tested_genes = []

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df)

    # log random number seed choice if provided
    if opts['seed'] is not None:
        logger.info('Pseudo Random Number Generator Seed: {0}'.format(opts['seed']))

    # perform permutation test
    bed_dict = utils.read_bed(opts['bed'], non_tested_genes)
    permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)

    # Perform BH p-value adjustment and tidy up data for output
    if opts['kind'] == 'oncogene':
        permutation_df = handle_oncogene_results(permutation_result, non_tested_genes)
    elif opts['kind'] == 'tsg':
        permutation_df = handle_tsg_results(permutation_result)
    elif opts['kind'] == 'effect':
        permutation_df = handle_effect_results(permutation_result)

    # save output
    if opts['output']:
        permutation_df.to_csv(opts['output'], sep='\t', index=False)

    return permutation_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
