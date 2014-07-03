# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))

# import extract_genes module
import extract_gene_seq as eg

def test_rev_comp():
    seq1 = 'CT'
    seq2 = 'AaCg'
    seq3 = 'aNnC'

    rc_seq1 = eg.rev_comp(seq1)
    assert rc_seq1 == 'AG'
    rc_seq2 = eg.rev_comp(seq2)
    assert rc_seq2 == 'cGtT'
    rc_seq3 = eg.rev_comp(seq3)
    assert rc_seq3 == 'GnNt'


def test_main():
    opts = {'input': os.path.join(file_dir, '../data/hg19.fa'),
            'output': os.path.join(file_dir, 'output/example_genes.fa'),
            'bed': os.path.join(file_dir, 'data/example.bed')}
    eg.main(opts)
