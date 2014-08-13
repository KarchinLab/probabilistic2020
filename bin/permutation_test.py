#!/usr/bin/env python
import sys
try:
    # fix problems with pythons terrible import system
    import os
    file_dir = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(os.path.join(file_dir, '../permutation2020/python'))
    sys.path.append(os.path.join(file_dir, '../permutation2020/cython'))

    # normal imports
    import utils
    from gene_sequence import GeneSequence
    from sequence_context import SequenceContext
    import cutils
except:
    import permutation2020.python.utils as utils
    from permutation2020.python.gene_sequence import GeneSequence
    from permutation2020.python.sequence_context import SequenceContext
    import permutation2020.cython.cutils as cutils

import argparse
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool
import logging
import datetime
import itertools as it
from functools import wraps

logger = logging.getLogger(__name__)  # module logger

def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def log_error_decorator(f):
    """Writes exception to log file if occured in decorated function.

    This decorator wrapper is needed for multiprocess logging since otherwise
    the python multiprocessing module will obscure the actual line of the error.
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        try:
            result = f(*args, **kwds)
            return result
        except KeyboardInterrupt:
            logger.info('Ctrl-C stopped a process.')
        except Exception, e:
            logger.exception(e)
            raise
    return wrapper


def keyboard_exit_wrapper(func):
    def wrap(self, timeout=None):
        # Note: the timeout of 1 googol seconds introduces a rather subtle
        # bug for Python scripts intended to run many times the age of the universe.
        return func(self, timeout=timeout if timeout is not None else 1e100)
    return wrap


def deleterious_permutation(context_counts,
                            context_to_mut,
                            seq_context,
                            gene_seq,
                            num_permutations=10000):
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
        tmp_mut_info = utils.get_aa_mut_info(row,
                                             somatic_base,
                                             gene_seq)

        # calc deleterious mutation info
        tmp_del_count = cutils.calc_deleterious_info(tmp_mut_info['Reference AA'],
                                                     tmp_mut_info['Somatic AA'])
        del_count_list.append(tmp_del_count)
    return del_count_list


def position_permutation(context_counts,
                         context_to_mut,
                         seq_context,
                         gene_seq,
                         num_permutations=10000,
                         kde_bandwidth=None):
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
    kde_bandwidth : int, default: None
        ?possibly deprecated parameter

    Returns
    -------
    num_recur_list : list
        list of recurrent mutation counts under the null
    entropy_list : list
        list of position entropy values under the null
    entropy_list : list
        list of position entropy values after KDE smoothing
        under the null.
    bw_list : list
        list of cross-validated KDE bandwidth values under the null.
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
    num_recur_list, entropy_list, kde_entropy_list, bw_list = [], [], [], []
    for row in tmp_mut_pos:
        # get info about mutations
        tmp_mut_info = utils.get_aa_mut_info(row,
                                             somatic_base,
                                             gene_seq)

        # calculate position info
        tmp_recur_ct, tmp_entropy, tmp_kde_ent, tmp_bw = cutils.calc_pos_info(tmp_mut_info['Codon Pos'],
                                                                              tmp_mut_info['Reference AA'],
                                                                              tmp_mut_info['Somatic AA'],
                                                                              kde_bandwidth)
        num_recur_list.append(tmp_recur_ct)
        entropy_list.append(tmp_entropy)
        kde_entropy_list.append(tmp_kde_ent)
        bw_list.append(tmp_bw)

    return num_recur_list, entropy_list, kde_entropy_list, bw_list


def read_bed(file_path, filtered_genes):
    """Reads BED file and populates a dictionary separating genes
    by chromosome.

    Parameters
    ----------
    file_path : str
        path to BED file
    filtered_genes: list
        list of gene names to not use

    Returns
    -------
    bed_dict: dict
        dictionary mapping chromosome keys to a list of BED lines
    """
    # read in entire bed file into a dict with keys as chromsomes
    bed_dict = {}
    for bed_row in utils.bed_generator(file_path):
        if bed_row.gene_name not in filtered_genes:
            bed_dict.setdefault(bed_row.chrom, [])
            bed_dict[bed_row.chrom].append(bed_row)
    return bed_dict


def calc_deleterious_p_value(mut_info,
                             unmapped_mut_info,
                             sc,
                             gs,
                             bed,
                             num_permutations,
                             del_threshold):
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
        aa_mut_info = utils.get_aa_mut_info(mut_info['Coding Position'],
                                            mut_info['Tumor_Allele'].tolist(),
                                            gs)
        ref_aa = aa_mut_info['Reference AA'] + unmapped_mut_info['Reference AA']
        somatic_aa = aa_mut_info['Somatic AA'] + unmapped_mut_info['Somatic AA']
        num_del = cutils.calc_deleterious_info(ref_aa, somatic_aa)

        # skip permutation test if number of deleterious mutations is not at
        # least meet some user-specified threshold
        if num_del >= del_threshold:
            # perform permutations
            null_del_list = deleterious_permutation(context_cts,
                                                    context_to_mutations,
                                                    sc,  # sequence context obj
                                                    gs,
                                                    num_permutations)  # gene sequence obj

            # calculate p-value
            del_num_nulls = sum([1 for d in null_del_list
                                 if d >= num_del])
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
                          num_permutations):
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
        permutation_result = position_permutation(context_cts,
                                                  context_to_mutations,
                                                  sc,  # sequence context obj
                                                  gs,
                                                  num_permutations)  # gene sequence obj
        num_recur_list, pos_entropy_list, kde_entropy_list, bw_list = permutation_result  # unpack results

        # get recurrent info for actual mutations
        aa_mut_info = utils.get_aa_mut_info(mut_info['Coding Position'],
                                            mut_info['Tumor_Allele'].tolist(),
                                            gs)
        codon_pos = aa_mut_info['Codon Pos'] + unmapped_mut_info['Codon Pos']
        ref_aa = aa_mut_info['Reference AA'] + unmapped_mut_info['Reference AA']
        somatic_aa = aa_mut_info['Somatic AA'] + unmapped_mut_info['Somatic AA']
        num_recurrent, pos_ent, kde_ent, opt_bw = cutils.calc_pos_info(codon_pos,
                                                                       ref_aa,
                                                                       somatic_aa,
                                                                       None)

        # calculate permutation p-value
        recur_num_nulls = sum([1 for null_recur in num_recur_list
                               if null_recur >= num_recurrent])
        entropy_num_nulls = sum([1 for null_ent in pos_entropy_list
                                 if null_ent <= pos_ent])
        kde_entropy_num_nulls = sum([1 for null_ent in kde_entropy_list
                                     if null_ent <= kde_ent])
        kde_bw_num_nulls = sum([1 for null_bw in bw_list
                                if null_bw <= opt_bw])
        recur_p_value = recur_num_nulls / float(num_permutations)
        ent_p_value = entropy_num_nulls / float(num_permutations)
        kde_ent_p_value = kde_entropy_num_nulls / float(num_permutations)
        kde_bw_p_value = kde_bw_num_nulls / float(num_permutations)
    else:
        num_recurrent = 0
        pos_ent = 0
        kde_ent = 0
        opt_bw = 0
        recur_p_value = 1.0
        ent_p_value = 1.0
        kde_ent_p_value = 1.0
        kde_bw_p_value = 1.0
    result = [bed.gene_name, num_recurrent, pos_ent, kde_ent, opt_bw,
              recur_p_value, ent_p_value, kde_ent_p_value, kde_bw_p_value]
    return result


def recover_unmapped_mut_info(mut_info, bed, sc, opts):
    # retreive info based on annotated protein effects and genomic coordinates
    if opts['use_unmapped'] and opts['genome']:
        genome_fa = pysam.Fastafile(opts['genome'])
        # try to still use mutations that are not on the reference transcript
        tmp_mut_info = mut_info[mut_info['Coding Position'].isnull()]
        unmapped_mut_info = utils.get_unmapped_aa_mut_info(tmp_mut_info,
                                                           genome_fa,
                                                           bed.strand,
                                                           bed.chrom,
                                                           opts['context'])
        genome_fa.close()

        # filter out cases where the nucleotide context does not exist
        # on the reference transcript
        bad_contexts = [i for i in range(len(unmapped_mut_info['Context']))
                        if not sc.is_valid_context(unmapped_mut_info['Context'][i])]
        for key in unmapped_mut_info:
            unmapped_mut_info[key] = utils.filter_list(unmapped_mut_info[key],
                                                       bad_contexts)
    else:
        unmapped_mut_info = {'Context': [], 'Reference AA': [], 'Codon Pos': [],
                             'Somatic AA': [], 'Tumor_Allele': []}
    return unmapped_mut_info


@log_error_decorator
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
        sc = SequenceContext(gs)

        # count total mutations in gene
        total_mut = len(mut_info)

        # fix nucleotide letter if gene is on - strand
        if bed.strand == '-':
            mut_info['Tumor_Allele'].map(lambda x: utils.rev_comp(x))

        # get coding positions, mutations unmapped to the reference tx will have
        # NA for a coding position
        pos_list = []
        for ix, row in mut_info.iterrows():
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], row['Start_Position'])
            pos_list.append(coding_pos)
        mut_info['Coding Position'] = pos_list

        # recover mutations that could not be mapped to the reference transcript
        # for a gene before being dropped (next step)
        unmapped_mut_info = recover_unmapped_mut_info(mut_info, bed, sc, opts)

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
                                               gs, bed, num_permutations)
            result.append(tmp_result + [total_mut, unmapped_muts])
        else:
            # calculate results for deleterious mutation permutation test
            tmp_result = calc_deleterious_p_value(mut_info, unmapped_mut_info,
                                                  sc, gs, bed, num_permutations,
                                                  opts['deleterious'])
            result.append(tmp_result + [total_mut, unmapped_muts])

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def multiprocess_permutation(bed_dict, mut_df, opts):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys())
    num_processes = opts['processes']
    result_list = []
    for i in range(0, len(chroms), num_processes):
        pool = Pool(processes=num_processes)
        tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
        info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                        for tmp_ix in range(i, i+tmp_num_proc))
        process_results = pool.imap(singleprocess_permutation, info_repeat)
        process_results.next = keyboard_exit_wrapper(process_results.next)
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
    col_order  = ['gene', 'Total Mutations', 'Unmapped to Ref Tx',
                  'num deleterious', 'deleterious p-value',
                  'deleterious BH q-value']
    return permutation_df[col_order]


def handle_oncogene_results(permutation_result, non_tested_genes):
    permutation_df = pd.DataFrame(sorted(permutation_result, key=lambda x: x[5]),
                                  columns=['gene', 'num recurrent', 'position entropy',
                                           'kde position entropy', 'kde bandwidth', 'recurrent p-value',
                                           'entropy p-value', 'kde entropy p-value', 'kde bandwidth p-value',
                                           'Total Mutations', "Unmapped to Ref Tx"])

    # get benjamani hochberg adjusted p-values
    permutation_df['recurrent BH q-value'] = utils.bh_fdr(permutation_df['recurrent p-value'])
    permutation_df['entropy BH q-value'] = utils.bh_fdr(permutation_df['entropy p-value'])
    permutation_df['kde entropy BH q-value'] = utils.bh_fdr(permutation_df['kde entropy p-value'])
    permutation_df['kde bandwidth BH q-value'] = utils.bh_fdr(permutation_df['kde bandwidth p-value'])

    # include non-tested genes in the result
    no_test_df = pd.DataFrame(index=range(len(non_tested_genes)))
    no_test_df['Performed Recurrency Test'] = 0
    no_test_df['gene'] = non_tested_genes
    permutation_df = pd.concat([permutation_df, no_test_df])
    permutation_df['Performed Recurrency Test'] = permutation_df['Performed Recurrency Test'].fillna(1).astype(int)

    # order output
    permutation_df['num recurrent'] = permutation_df['num recurrent'].fillna(-1).astype(int)  # fix dtype isssue
    col_order = ['gene', 'Total Mutations', 'Unmapped to Ref Tx', 'num recurrent', 'position entropy', 'kde position entropy', 'kde bandwidth',
                 'recurrent p-value', 'recurrent BH q-value', 'entropy p-value', 'entropy BH q-value',
                 'kde entropy p-value', 'kde entropy BH q-value', 'kde bandwidth p-value',
                 'kde bandwidth BH q-value', 'Performed Recurrency Test']
    return permutation_df[col_order]


def _fix_mutation_df(mutation_df):
    allowed_types = ['Missense_Mutation', 'Silent', 'Nonsense_Mutation', 'Splice_Site']
    mutation_df = mutation_df[mutation_df.Variant_Classification.isin(allowed_types)]  # only keep SNV
    valid_nuc_flag = (mutation_df['Reference_Allele'].apply(utils.is_valid_nuc) & \
                      mutation_df['Tumor_Allele'].apply(utils.is_valid_nuc))
    mutation_df = mutation_df[valid_nuc_flag]  # filter bad lines
    mutation_df['Start_Position'] = mutation_df['Start_Position'] - 1
    mutation_df = mutation_df[mutation_df['Tumor_Allele'].apply(lambda x: len(x)==1)]
    mutation_df = mutation_df[mutation_df['Reference_Allele'].apply(lambda x: len(x)==1)]
    return mutation_df


def _get_high_tsg_score(mutation_df, tsg_score_thresh):
    mutation_df['indicator'] = 1
    table = pd.pivot_table(mutation_df,
                           values='indicator',
                           cols='Variant_Classification',
                           rows='Gene',
                           aggfunc=np.sum)
    mut_type_frac = table.div(table.sum(axis=1).astype(float), axis=0).fillna(0.0)
    for c in ['Nonsense_Mutation', 'Frame_Shift_Indel', 'Splice_Site', 'Nonstop_Mutation']:
        if c not in mut_type_frac.columns:
            mut_type_frac[c] = 0.0  # make sure columns are defined
    tsg_score = mut_type_frac['Nonsense_Mutation'] + mut_type_frac['Frame_Shift_Indel'] + \
                mut_type_frac['Splice_Site'] + mut_type_frac['Nonstop_Mutation']
    non_tested_genes = set(tsg_score[tsg_score>=tsg_score_thresh].index.tolist())
    return non_tested_genes


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
    help_str = 'Number of processes to use (more==Faster, default: 1).'
    parser.add_argument('-p', '--processes',
                        type=int, default=1,
                        help=help_str)
    help_str = 'Number of permutations for null model'
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
    help_str = ('Perform recurrent mutation permutation test if gene has '
                'atleast a user specified number of recurrent mutations (default: 2)')
    parser.add_argument('-r', '--recurrent',
                        type=int, default=2,
                        help=help_str)
    help_str = ('Perform tsg permutation test if gene has '
                'at least a user specified number of deleterious mutations (default: 5)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=5,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. (Default: .10)')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=.10,
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
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging

    opts = vars(args)
    if opts['use_unmapped'] and not opts['genome']:
        print('You must specify a genome fasta with -g if you set the '
              '--use-unmapped flag to true.')
        sys.exit(1)
    return opts


def main(opts):
    # hack to index the FASTA file
    gene_fa = pysam.Fastafile(opts['input'])
    gene_fa.close()

    # Get Mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # specify genes to skip
    if opts['kind'] == 'oncogene':
        # find genes with tsg score above threshold to filter out for oncogene
        # permutation test
        non_tested_genes = _get_high_tsg_score(mut_df, opts['tsg_score'])
    else:
        # don't filter out genes for tsg permutation test
        non_tested_genes = []

    # select valid single nucleotide variants only
    mut_df = _fix_mutation_df(mut_df)

    # perform permutation test
    bed_dict = read_bed(opts['bed'], non_tested_genes)
    permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)

    # Perform BH p-value adjustment and tidy up data for output
    if opts['kind'] == 'oncogene':
        permutation_df = handle_oncogene_results(permutation_result, non_tested_genes)
    else:
        permutation_df = handle_tsg_results(permutation_result)

    # save output
    permutation_df.to_csv(opts['output'], sep='\t', index=False)

    return permutation_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
