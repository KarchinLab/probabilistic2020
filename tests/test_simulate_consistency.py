# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))

import simulate_consistency as sc


def test_100genes_main():
    # check oncogene simulation
    opts = {'input': os.path.join(file_dir, 'data/100genes.fa'),
            'bed': os.path.join(file_dir, 'data/100genes.bed'),
            'mutations': os.path.join(file_dir, 'data/100genes_mutations.txt'),
            'output': os.path.join(file_dir, 'output/100genes_position_sim_consis_chasm_output.txt'),
            'context': 1.5,
            'tsg_score': .1,
            'recurrent': 3,
            'fraction': .02,
            'recurrent_pseudo_count': 0,
            'deleterious_pseudo_count': 0,
            'iterations': 5,
            'use_unmapped': False,
            'processes': 0,
            'num_permutations': 1000,
            'kind': 'oncogene'}
    result = sc.main(opts)

    # check tsg simulation
    opts['deleterious'] = 1
    opts['kind'] = 'tsg'
    opts['output'] = os.path.join(file_dir, 'output/100genes_deleterious_sim_consis_chasm_output.txt')
    result = sc.main(opts)
