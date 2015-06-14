# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))

import permutation_test as pt
import numpy as np


def test_tp53_main():
    opts = {'input': os.path.join(file_dir, 'data/tp53.fa'),
            'bed': os.path.join(file_dir, 'data/tp53.bed'),
            'mutations': os.path.join(file_dir, 'data/tp53_mutations.txt'),
            'output': os.path.join(file_dir, 'output/tp53_output.txt'),
            'context': 1,
            'use_unmapped': False,
            'deleterious': 5,
            'processes': 0,
            'num_permutations': 10000,
            'deleterious_pseudo_count': 0,
            'seed': None,
            'kind': 'tsg'}
    # single nucleotide context
    result = pt.main(opts)
    assert result.ix[0, 'inactivating p-value'] < 0.001, 'TP53 should have a very low p-value ({0}>.001)'.format(result[0][2])

    # di-nucleotide case
    opts['context'] = 2
    result = pt.main(opts)
    assert result.ix[0, 'inactivating p-value'] < 0.001, 'TP53 should have a very low p-value ({0}>.001)'.format(result[0][2])

    # no context case
    opts['context'] = 0
    result = pt.main(opts)
    assert result.ix[0, 'inactivating p-value'] < 0.001, 'TP53 should have a very low p-value ({0}>.001)'.format(result[0][2])


def test_100genes_main():
    opts = {'input': os.path.join(file_dir, 'data/100genes.fa'),
            'bed': os.path.join(file_dir, 'data/100genes.bed'),
            'mutations': os.path.join(file_dir, 'data/100genes_mutations.txt'),
            'output': os.path.join(file_dir, 'output/100genes_deleterious_single_nuc_output.txt'),
            'context': 1,
            'use_unmapped': False,
            'deleterious': 5,
            'processes': 0,
            'num_permutations': 1000,
            'deleterious_pseudo_count': 0,
            'seed': None,
            'kind': 'tsg'}
    # single nucleotide context
    result = pt.main(opts)
    num_del_sig = np.sum(result['inactivating BH q-value'] < .1)
    assert num_del_sig < 7, 'Few of the 100 test genes should not be significant ({0})'.format(num_del_sig)

    # no context case
    opts['context'] = 0
    opts['output'] = os.path.join(file_dir, 'output/100genes_deleterious_no_context_output.txt')
    result = pt.main(opts)
    num_del_sig = np.sum(result['inactivating BH q-value'] < .1)
    assert num_del_sig < 9, 'Few of the 100 test genes should not be significant ({0})'.format(num_del_sig)

    # di-nucleotide context
    opts['context'] = 2
    opts['output'] = os.path.join(file_dir, 'output/100genes_deleterious_dinuc_output.txt')
    result = pt.main(opts)
    num_del_sig = np.sum(result['inactivating BH q-value'] < .1)
    assert num_del_sig < 7, 'Few of the 100 test genes should not be significant ({0})'.format(num_del_sig)
