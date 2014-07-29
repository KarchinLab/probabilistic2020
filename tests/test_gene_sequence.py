# fix problems with pythons terrible import system
import os
import sys
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../bin/'))
sys.path.append(os.path.join(file_dir, '../src/python/'))
sys.path.append(os.path.join(file_dir, '../src/cython/'))

# useful imports
from gene_sequence import GeneSequence
import pysam
import utils
import cutils

# read in fake sequence
fake_fasta = os.path.join(file_dir, 'data/fake_sequence.fa')
fake_bed = os.path.join(file_dir, 'data/fake_gene.bed')
gene_fa = pysam.Fastafile(fake_fasta)
with open(fake_bed) as handle:
    bed = utils.BedLine(handle.readline().strip())

def test_simple_constructor():
    gs = GeneSequence(gene_fa, nuc_context=1)
    gs.set_gene(bed)
    assert gs.exon_seq == 'ACATGAATGATAGATCCGAAA', 'Sequence is not correct'

    # this should update the sequence correctly
    fake_germline = ['A', 'C', 'N', 'G', 'T']
    fake_pos = [1, 0, 20, 7, 15]
    gs.add_germline_variants(fake_germline, fake_pos)
    assert gs.exon_seq == 'CAATGAAGGATAGATTCGAAN'


def test_pos_to_codon():
    gs = GeneSequence(gene_fa, nuc_context=1)
    gs.set_gene(bed)

    pos_list = [1, 12, 17]
    results = []
    for pos in pos_list:
        codon_info = cutils.pos_to_codon(gs.exon_seq, pos)
        results.append(codon_info)
    true_results = [('ACA', 0, 1), ('GAT', 4, 0), ('CCG', 5, 2)]
    assert results == true_results, 'Codon information is incorrect'
