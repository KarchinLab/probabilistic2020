# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import prob2020.console.count_frameshifts as cf

def test_tp53_main():
    opts = {'mutations': os.path.join(file_dir, 'data/tp53_fs_mutations.txt'),
            'bins': 1,
            'sample_number': 294,
            'use_unmapped': False,
            'bed': os.path.join(file_dir, 'data/tp53.bed'),
            'output': os.path.join(file_dir, 'output/tp53_fs_counts.txt'),
           }
    fs_df = cf.main(opts)

    assert fs_df.loc['TP53', '1'] == 295, "Number of frameshifts should equal to 295"
