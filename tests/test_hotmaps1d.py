# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '..'))

import prob2020.console.randomization_test as rt
import numpy as np

def test_ctnnb1_hotmaps_main():
    opts = {'input': os.path.join(file_dir, 'data/CTNNB1.fa'),
            'bed': os.path.join(file_dir, 'data/CTNNB1.bed'),
            'mutations': os.path.join(file_dir, 'data/CTNNB1_mutations.txt'),
            'output': os.path.join(file_dir, 'output/CTNNB1_output_hotmaps.txt'),
            'context': 1.5,
            'use_unmapped': False,
            'processes': 0,
            'num_iterations': 1000,
            'stop_criteria': 100,
            'unique': 0,
            'seed': None,
            'window': '3',
            'report_index': True,
            'null_distr_dir': os.path.join(file_dir, 'output/hotmaps1d_null'),
            'kind': 'hotmaps1d'}
    # single nucleotide context
    result = rt.main(opts)

    # di-nucleotide case
    opts['window'] = '6'
    result = rt.main(opts)

    # no context case
    opts['window'] = '9'
    result = rt.main(opts)


def test_100genes_main():
    opts = {'input': os.path.join(file_dir, 'data/100genes.fa'),
            'bed': os.path.join(file_dir, 'data/100genes.bed'),
            'mutations': os.path.join(file_dir, 'data/100genes_mutations.txt'),
            'output': os.path.join(file_dir, 'output/100genes_hotmaps_single_nuc_output.txt'),
            'context': 1,
            'use_unmapped': False,
            'processes': 0,
            'num_iterations': 1000,
            'stop_criteria': 100,
            'unique': False,
            'seed': None,
            'window': '3',
            'report_index': False,
            'null_distr_dir': os.path.join(file_dir, 'output/hotmaps1d_null'),
            'kind': 'hotmaps1d'}
    # single nucleotide context
    result = rt.main(opts)
    num_sig = np.sum(result['q-value'] < .01)
    assert num_sig < 9, 'Few of the 100 test genes should not be significant ({0})'.format(num_sig)


if __name__ == '__main__':
    test_100genes_main()
