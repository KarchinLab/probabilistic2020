# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import prob2020.console.annotate as sm

def test_sim_summary():
    opts = {'input': os.path.join(file_dir, 'data/sim_summary.fa'),
            'mutations': os.path.join(file_dir, 'data/sim_summary_mutations.txt'),
            'bed': os.path.join(file_dir, 'data/sim_summary.bed'),
            'processes': 0,
            'num_iterations': 1,
            'context': 1.5,
            'summary': False,
            'maf': True,
            'unique': True,
            'use_unmapped': False,
            'genome': '',
            'score_dir': None,  # just skip using scores
            'fraction': .02,
            'recurrent': 3,
            'seed': 101,
            'output': os.path.join(file_dir, 'output/sim_summary_maf.txt')
            }
    sm.main(opts)

    opts['maf'] = False
    opts['summary'] = True
    opts['output'] = os.path.join(file_dir, 'output/sim_summary.txt')
    sm.main(opts)


if __name__ == '__main__':
    test_sim_summary()
