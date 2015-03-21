# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))
sys.path.append(os.path.join(file_dir, '../prob2020/python/'))
sys.path.append(os.path.join(file_dir, '../prob2020/cython/'))

# useful imports
import permutation_test as pt
import cutils
import utils
import numpy as np
import pandas as pd


def test_kde_entropy():
    import pysam
    from gene_sequence import GeneSequence

    # read fasta
    ctnnb1_fasta = os.path.join(file_dir, 'data/CTNNB1.fa')
    gene_fa = pysam.Fastafile(ctnnb1_fasta)
    gs = GeneSequence(gene_fa, nuc_context=1)

    # read CTNNB1 bed file
    ctnnb1_bed = os.path.join(file_dir, 'data/CTNNB1.bed')
    bed_list = [b for b in utils.bed_generator(ctnnb1_bed)]
    gs.set_gene(bed_list[0])

    # specify mutation
    coding_pos = [0]
    somatic_base = ['C']

    # check mutation info
    aa_info = utils.get_aa_mut_info(coding_pos, somatic_base, gs)
    pos_array = np.array(aa_info['Codon Pos'], dtype=np.int)
    entropy, bandwidth = cutils.kde_entropy(pos_array)
    assert_msg = ('A single mutation should be 1.0 for entropy fraction '
                  '(entropy fraction: {0}, bandwidth: {1})'.format(entropy, bandwidth))
    assert entropy == 1.0, assert_msg

    opts = {'input': os.path.join(file_dir, 'data/CTNNB1.fa'),
            'bed': os.path.join(file_dir, 'data/CTNNB1.bed'),
            'mutations': os.path.join(file_dir, 'data/CTNNB1_mutations.txt'),
            'output': os.path.join(file_dir, 'output/CTNNB1_output.txt'),
            'context': 1,
            'tsg_score': .05,
            'processes': 1,
            'num_permutations': 10000,
            'kind': 'oncogene'}
    mut_df = pd.read_csv(opts['mutations'], sep='\t')

    # CTNNB1 should have few deleterious mutations, so check it
    non_tested_genes = pt._get_high_tsg_score(mut_df, opts['tsg_score'])
    assert len(non_tested_genes)==0, 'CTNNB1 should have very few deleterious mutations'

    # perform permutation test only on CTNNB1
    mut_df = pt._fix_mutation_df(mut_df)
    info = (bed_list, mut_df, opts)  # singleprocess_permutation takes a single tuple argument
    perm_result = pt.singleprocess_permutation(info)
    raise ValueError(perm_result)
