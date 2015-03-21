# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))
sys.path.append(os.path.join(file_dir, '..'))

# useful imports
from prob2020.python.gene_sequence import GeneSequence
from prob2020.python.sequence_context import SequenceContext
import prob2020.python.utils as utils
import pysam

# set up global variables
fake_fasta = os.path.join(file_dir, 'data/fake_sequence.fa')
fake_bed = os.path.join(file_dir, 'data/fake_gene.bed')
gene_fa = pysam.Fastafile(fake_fasta)
with open(fake_bed) as handle:
    bed = utils.BedLine(handle.readline().strip())


def test_single_context_constructor():
    # single nuc context
    gs = GeneSequence(gene_fa, nuc_context=1)
    gs.set_gene(bed)
    sc = SequenceContext(gs)
    true_counts = {'A': 10,
                   'T': 4,
                   'G': 4,
                   'C': 3}
    true_ctxt2pos = {'A': [0, 2, 5, 6, 9, 11, 13, 18, 19, 20],
                     'C': [1, 15, 16],
                     'G': [4, 8, 12, 17],
                     'T': [3, 7, 10, 14]}
    _check_true_counts(sc, true_counts)
    _check_true_context_pos(sc, true_ctxt2pos)


def test_dinuc_context_constructor():
    # dinucleotide context
    gs = GeneSequence(gene_fa, nuc_context=2)
    gs.set_gene(bed)
    sc = SequenceContext(gs)


def test_chasm_context_constructor():
    # chasm context, mixture between single and di context
    gs = GeneSequence(gene_fa, nuc_context=1.5)
    gs.set_gene(bed)
    sc = SequenceContext(gs)
    true_counts = {'A': 10,
                   'C': 1,
                   'C*pG': 1,
                   'CpG*': 1,
                   'G*pA': 3,
                   'T': 4,
                   'TpC*': 1}
    true_ctxt2pos = {'A': [2, 5, 6, 9, 11, 13, 18, 19, 0, 20],
                     'C': [1],
                     'C*pG': [16],
                     'CpG*': [17],
                     'G*pA': [4, 8, 12],
                     'T': [3, 7, 10, 14],
                     'TpC*': [15]}
    _check_true_counts(sc, true_counts)
    _check_true_context_pos(sc, true_ctxt2pos)


def test_no_context_constructor():
    # no context
    gs = GeneSequence(gene_fa, nuc_context=0)
    gs.set_gene(bed)
    sc = SequenceContext(gs)
    true_counts = {'None': 21}
    true_ctxt2pos = {'None': range(21)}
    _check_true_counts(sc, true_counts)
    _check_true_context_pos(sc, true_ctxt2pos)


def _check_true_counts(seq_context, true_counts):
    for letter in true_counts:
        true_ct = true_counts[letter]
        ct = len(seq_context.context2pos[letter])
        assert true_ct == ct, 'Context count mismatch'


def _check_true_context_pos(seq_context, true_context_pos):
    for letter in true_context_pos:
        assert_msg = 'Context positions don\'t match ({0}: {1} != {2})'.format(letter,
                                                                               true_context_pos[letter],
                                                                               seq_context.context2pos[letter])
        assert true_context_pos[letter] == seq_context.context2pos[letter], assert_msg
