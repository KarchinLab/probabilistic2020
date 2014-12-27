# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))

import simulate_summary as sm

def test_sim_summary():
    opts = {'input': os.path.join(file_dir, 'data/sim_summary.fa'),
            'mutations': os.path.join(file_dir, 'data/sim_summary_mutations.txt'),
            'bed': os.path.join(file_dir, 'data/sim_summary.bed'),
            'processes': 0,
            'num_permutations': 1,
            'context': 1.5,
            'summary': False,
            'maf': True,
            'use_unmapped': False,
            'genome': '',
            'output': os.path.join(file_dir, 'output/sim_summary_maf.txt')
            }
    sm.main(opts)

    opts['maf'] = False
    opts['summary'] = True
    opts['output'] = os.path.join(file_dir, 'output/sim_summary.txt')
    sm.main(opts)


if __name__ == '__main__':
    test_sim_summary()
