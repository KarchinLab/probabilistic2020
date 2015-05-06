import numpy as np
import prob2020.python.mutation_context as mc
import prob2020.python.utils as utils
from ..cython import cutils

def annotate_maf(coding_pos, somatic_base, gene_seq):
    # make sure numpy array
    coding_pos = np.array(coding_pos)

    # info about gene
    gene_name = gene_seq.bed.gene_name
    strand = gene_seq.bed.strand
    chrom = gene_seq.bed.chrom
    gene_seq.bed.init_genome_coordinates()  # map seq pos to genome

    # determine result of random positions
    maf_list = []

    # get genome coordinate
    pos2genome = np.vectorize(lambda x: gene_seq.bed.seqpos2genome[x]+1)
    genome_coord = pos2genome(coding_pos)

    # get info about mutations
    tmp_mut_info = mc.get_aa_mut_info(coding_pos,
                                      somatic_base,
                                      gene_seq)

    # get string describing variant
    var_class = cutils.get_variant_classification(tmp_mut_info['Reference AA'],
                                                  tmp_mut_info['Somatic AA'],
                                                  tmp_mut_info['Codon Pos'])

    # prepare output
    for k, mysomatic_base in enumerate(somatic_base):
        ######
        # Note: positions are converted to 1-based positions
        # for reporting DNA/Protein change, but internally
        # they are represented as 0-based
        ######
        # format DNA change
        ref_nuc = tmp_mut_info['Reference Nuc'][k]
        nuc_pos = coding_pos[k]
        dna_change = 'c.{0}{1}>{2}'.format(ref_nuc, nuc_pos+1, mysomatic_base)

        # format protein change
        ref_aa = tmp_mut_info['Reference AA'][k]
        somatic_aa = tmp_mut_info['Somatic AA'][k]
        codon_pos = tmp_mut_info['Codon Pos'][k]
        codon_pos_1_based = (codon_pos + 1) if codon_pos is not None else None
        protein_change = 'p.{0}{1}{2}'.format(ref_aa, codon_pos_1_based, somatic_aa)

        # reverse complement if on negative strand
        if strand == '-':
            ref_nuc = utils.rev_comp(ref_nuc)
            mysomatic_base = utils.rev_comp(mysomatic_base)

        # append results
        maf_line = [gene_name, strand, chrom, genome_coord[k], genome_coord[k],
                    ref_nuc, mysomatic_base, dna_change,
                    protein_change, var_class[k]]
        maf_list.append(maf_line)

    return maf_list
