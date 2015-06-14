# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))

import simulate_performance as sp


def test_100genes_main():
    # check oncogene simulation
    oncogene_path = 'output/100genes_position_sim_perf_chasm_output.txt'
    opts = {'input': os.path.join(file_dir, 'data/100genes.fa'),
            'bed': os.path.join(file_dir, 'data/100genes.bed'),
            'mutations': os.path.join(file_dir, 'data/100genes_mutations.txt'),
            'non_coding_background': os.path.join(file_dir, 'data/non_coding_fs.background.txt'),
            'bins': 5,
            'sample_number': 9000,
            'text_output': os.path.join(file_dir, oncogene_path),
            'plot_output': os.path.join(file_dir, 'output/'),
            'output': '',
            'context': 1.5,
            'tsg_score': .1,
            'recurrent': 3,
            'fraction': .02,
            'recurrent_pseudo_count': 0,
            'deleterious_pseudo_count': 0,
            'iterations': 5,
            'use_unmapped': False,
            'processes': 0,
            'num_permutations': 100,
            'kind': 'oncogene',
            'random_samples': True,
            'random_tumor_types': False,
            'bootstrap': False,
            'seed': None,
            'start_sample_rate': .1,
            'end_sample_rate': 1.0,
            'num_sample_rate': 3}
    result = sp.main(opts)

    # check tsg simulation
    opts['deleterious'] = 1
    opts['kind'] = 'tsg'
    tsg_path = 'output/100genes_deleterious_sim_perf_chasm_output.txt'
    opts['text_output'] = os.path.join(file_dir, tsg_path)
    result = sp.main(opts)

    # check effect simulation
    opts['kind'] = 'effect'
    effect_path = 'output/100genes_effect_sim_perf_chasm_output.txt'
    opts['text_output'] = os.path.join(file_dir, effect_path)
    result = sp.main(opts)
