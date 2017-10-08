#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

# package import
import prob2020.python.utils as utils
from prob2020.python.gene_sequence import GeneSequence
import prob2020.cython.cutils as cutils
import prob2020.python.mutation_context as mc
import prob2020.python.permutation as pm
import prob2020.python.indel as indel
import prob2020.python.annotate as anot
import prob2020.python.mymath as math

# external imports
import numpy as np
import pandas as pd
import pysam
import csv
from multiprocessing import Pool
import argparse
import logging
import copy
import itertools as it

logger = logging.getLogger(__name__)  # module logger

def multiprocess_permutation(bed_dict, mut_df, opts, indel_df=None):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys(), key=lambda x: len(bed_dict[x]), reverse=True)
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    #file_handle = open(opts['output'], 'w')
    file_handle = opts['handle']
    mywriter = csv.writer(file_handle, delimiter='\t', lineterminator='\n')
    if opts['maf'] and opts['num_iterations']:
        header = ['Gene', 'strand', 'Chromosome', 'Start_Position',
                  'End_Position', 'Reference_Allele', 'Tumor_Allele',
                  'Context', 'DNA_Change', 'Protein_Change', 'Variant_Classification']
    elif opts['maf']:
        header = ['Gene', 'strand', 'Chromosome', 'Start_Position',
                  'End_Position', 'Reference_Allele', 'Tumor_Allele',
                  'DNA_Change', 'Protein_Change', 'Variant_Classification',
                  'Tumor_Sample', 'Tumor_Type']
    else:
        header = ['Gene', 'ID', 'gene length', 'non-silent snv', 'silent snv', 'nonsense', 'lost stop',
                  'splice site', 'lost start', 'missense', 'recurrent missense',
                  'normalized missense position entropy',]
        # add column header for scores, is user provided one
        if opts['score_dir']:
            header += ['Total Missense MGAEntropy', 'Total Missense VEST Score']
        # add indel columns
        header += ['frameshift indel', 'inframe indel', 'normalized mutation entropy']
    mywriter.writerow(header)
    num_iterations = opts['num_iterations']

    # simulate indel counts
    if opts['summary'] and num_iterations:
        fs_cts, inframe_cts, gene_names = indel.simulate_indel_counts(indel_df,
                                                                      bed_dict,
                                                                      num_iterations,
                                                                      opts['seed'])
        name2ix = {gene_names[z]: z for z in range(len(gene_names))}
    # just count observed indels
    elif opts['summary']:
        # get gene names
        gene_names = [mybed.gene_name
                      for chrom in bed_dict
                      for mybed in bed_dict[chrom]]
        name2ix = {gene_names[z]: z for z in range(len(gene_names))}

        # initiate count vectors
        inframe_cts = np.zeros((1, len(gene_names)))
        fs_cts = np.zeros((1, len(gene_names)))

        # populate observed counts
        indel_cts_dict = indel_df['Gene'].value_counts().to_dict()
        fs_cts_dict = indel_df[indel.is_frameshift_annotation(indel_df)]['Gene'].value_counts().to_dict()
        for mygene in indel_cts_dict:
            if mygene in name2ix:
                # gene should be found in BED file annotation
                ix = name2ix[mygene]
                fs_cts[0, ix] = 0 if mygene not in fs_cts_dict else fs_cts_dict[mygene]
                inframe_cts[0, ix] = indel_cts_dict[mygene] - fs_cts[0, ix]

    # simulate snvs
    obs_result = []
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                # iterate through each chromosome result
                for chrom_result in process_results:
                    # add columns for indels
                    if opts['summary']:
                        tmp_chrom_result = []
                        for gname, grp in it.groupby(chrom_result, lambda x: x[0]):
                            for l, row in enumerate(grp):
                                gene_ix = name2ix[gname]
                                fs_count = fs_cts[l, gene_ix]
                                inframe_count = inframe_cts[l, gene_ix]
                                missense_pos_ct = list(row.pop(-1).values())  # missense codon counts
                                silent_pos_ct = [1 for l in range(row[4])]
                                inactivating_ct = sum(row[5:9]) + fs_count
                                tmp_count_list = missense_pos_ct + silent_pos_ct + [inactivating_ct, inframe_count]
                                norm_ent = math.normalized_mutation_entropy(tmp_count_list)
                                tmp_chrom_result.append(row+[fs_count, inframe_count, norm_ent])
                        chrom_result = tmp_chrom_result

                    # write output to file
                    mywriter.writerows(chrom_result)
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            # perform simulation
            info = (bed_dict[chroms[i]], mut_df, opts)
            chrom_results = singleprocess_permutation(info)

            # add indel columns
            if opts['summary']:
                tmp_chrom_result = []
                for gname, grp in it.groupby(chrom_results, lambda x: x[0]):
                    for l, row in enumerate(grp):
                        gene_ix = name2ix[gname]
                        fs_count = fs_cts[l, gene_ix]
                        inframe_count = inframe_cts[l, gene_ix]
                        missense_pos_ct = list(row.pop(-1).values())  # missense codon counts
                        silent_pos_ct = [1 for l in range(row[4])]
                        inactivating_ct = sum(row[5:9]) + fs_count
                        tmp_count_list = missense_pos_ct + silent_pos_ct + [inactivating_ct, inframe_count]
                        norm_ent = math.normalized_mutation_entropy(tmp_count_list)
                        tmp_chrom_result.append(row+[fs_count, inframe_count, norm_ent])
                chrom_results = tmp_chrom_result

            # write to file
            mywriter.writerows(chrom_results)
    #file_handle.close()


@utils.log_error_decorator
def singleprocess_permutation(info):
    bed_list, mut_df, opts = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    num_iterations = opts['num_iterations']
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    # go through each gene to perform simulation
    result = []
    for bed in bed_list:
        # compute context counts and somatic bases for each context
        gene_tuple = mc.compute_mutation_context(bed, gs, mut_df, opts)
        context_cts, context_to_mutations, mutations_df, gs, sc = gene_tuple

        if context_to_mutations:
            ## get information about observed non-silent counts
            if opts['summary'] and not num_iterations:
                tmp_mut_info = mc.get_aa_mut_info(mutations_df['Coding Position'],
                                                  mutations_df['Tumor_Allele'].tolist(),
                                                  gs)
                # calc mutation info summarizing observed mutations
                tmp_result = cutils.calc_summary_info(tmp_mut_info['Reference AA'],
                                                      tmp_mut_info['Somatic AA'],
                                                      tmp_mut_info['Codon Pos'],
                                                      bed.gene_name,
                                                      opts['score_dir'],
                                                      min_frac=opts['fraction'],
                                                      min_recur=opts['recurrent'])
                tmp_result = [[bed.gene_name, 'NA', bed.cds_len] + tmp_result]
            ## Just record protein changes in MAF
            elif opts['maf'] and not num_iterations:
                # input code for just annotating genes mutations
                tmp_result = anot.annotate_maf(mutations_df['Coding Position'],
                                               mutations_df['Tumor_Allele'].tolist(),
                                               gs)
                # add tumor sample / tumor type info to output
                tmp_result = [line + [mutations_df['Tumor_Sample'].iloc[i],
                                      mutations_df['Tumor_Type'].iloc[i]]
                              for i, line in enumerate(tmp_result)]
            ## Do permutations
            elif opts['maf']:
                # if user specified MAF format then output all mutations in
                # MAF format
                tmp_result = pm.maf_permutation(context_cts,
                                                context_to_mutations,
                                                sc,
                                                gs,
                                                num_iterations,
                                                drop_silent=opts['drop_silent'])
            else:
                # Summarized results for feature for each simulation for each
                # gene
                tmp_result = pm.summary_permutation(context_cts,
                                                    context_to_mutations,
                                                    sc,  # sequence context obj
                                                    gs,  # gene sequence obj
                                                    opts['score_dir'],
                                                    num_iterations,
                                                    min_frac=opts['fraction'],
                                                    min_recur=opts['recurrent'],
                                                    drop_silent=opts['drop_silent'])
            result += tmp_result

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def parse_arguments():
    # make a parser
    info = 'Either simulates or summarizes observed mutation data.'
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
                        default='stdout',
                        help='Path to log file. (accepts "stdout")')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

    # program arguments
    help_str = 'gene FASTA file from extract_gene_seq script'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'DNA mutations file (MAF file)'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'BED file annotation of genes'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Directory containing pre-compute score information in '
                'for VEST and evolutionary conservation in pickle format (Default: None).')
    parser.add_argument('-s', '--score-dir',
                        type=str, default=None,
                        help=help_str)
    help_str = ('Number of processes to use. 0 indicates using a single '
                'process without using a multiprocessing pool '
                '(more means Faster, default: 0).')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Number of iterations for null model simulations. If zero is '
                'specified then output represents a result from actually observed mutations (provided by -m parameter), '
                'otherwise results will be from simulated mutations. (Default: 0).')
    parser.add_argument('-n', '--num-iterations',
                        type=int, default=0,
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
    parser_grouper = parser.add_mutually_exclusive_group(required=True)
    parser_grouper.add_argument('--summary',
                                action='store_true',
                                help='Flag for saving results as summarized '
                                'features used (Default: True).')
    parser_grouper.add_argument('--maf',
                                action='store_true',
                                help='Flag for saving results in MAF format '
                                '(Default: False).')
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
                'as recurrently mutated (Default: .02).')
    parser.add_argument('-f', '--fraction',
                        type=float, default=.02,
                        help=help_str)
    help_str = ('Only keep unique mutations for each tumor sample.'
                'Mutations reproted from heterogeneous sources may contain'
                ' duplicates, e.g. a tumor sample was sequenced twice.')
    parser.add_argument('--unique',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Drop silent mutations in simulations. Useful if provided mutations '
                "don't include silent mutations")
    parser.add_argument('--drop-silent',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Important option for gene panels. Restrict simulated indels to only occur within the same set of '
                "genes as provied in the mutation file. ")
    parser.add_argument('--restrict-genes',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Specify the seed for the pseudo random number generator. '
                'By default, the seed is randomly chosen based. The seed will '
                'be used for the monte carlo simulations (Default: 101).')
    parser.add_argument('-seed', '--seed',
                        type=int, default=101,
                        help=help_str)
    help_str = 'Output text file of results'
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
                        log_level=log_level,
                        verbose=args.verbose)  # start logging

    opts = vars(args)
    if opts['use_unmapped'] and not opts['genome']:
        print('You must specify a genome fasta with -g if you set the '
              '--use-unmapped flag to true.')
        sys.exit(1)

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts):
    # hack to index the FASTA file
    gene_fa = pysam.Fastafile(opts['input'])
    gene_fa.close()

    # Get Mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)

    # rename columns to fit my internal column names
    rename_dict = {
        'Hugo_Symbol': 'Gene',
        'Tumor_Sample_Barcode': 'Tumor_Sample',
        'Tumor_Seq_Allele2' : 'Tumor_Allele'
    }
    mut_df.rename(columns=rename_dict, inplace=True)

    # restrict to only observed genes if flag present
    restricted_genes = None
    if opts['restrict_genes']:
        restricted_genes = set(mut_df['Gene'].unique())

    # process indels
    indel_df = indel.keep_indels(mut_df)  # return indels only
    indel_df.loc[:, 'Start_Position'] = indel_df['Start_Position'] - 1  # convert to 0-based
    indel_df.loc[:, 'indel len'] = indel_df['indel len'] + 1
    logger.info('There were {0} indels identified.'.format(len(indel_df)))
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df, opts['unique'])

    # read in bed info
    bed_dict = utils.read_bed(opts['bed'], restricted_genes)

    # perform permutation
    opts['handle'] = open(opts['output'], 'w')
    multiprocess_permutation(bed_dict, mut_df, opts, indel_df)

    # save indels
    if opts['maf']:
        #with open(opts['output'], 'a') as handle:
        mywriter = csv.writer(opts['handle'], delimiter='\t', lineterminator='\n')
        for maf_lines in indel.simulate_indel_maf(indel_df, bed_dict,
                                                  opts['num_iterations'],
                                                  opts['seed']):
            mywriter.writerows(maf_lines)
    opts['handle'].close()


def cli_main():
    opts = parse_arguments()
    main(opts)


if __name__ == "__main__":
    cli_main()
