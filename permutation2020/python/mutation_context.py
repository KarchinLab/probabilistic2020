import permutation2020.python.utils as utils
import permutation2020.python.sequence_context
from amino_acid import AminoAcid
from ..cython import cutils
import numpy as np
import pandas as pd
import pysam
import itertools as it


def get_all_context_names(context_num):
    if context_num == 0:
        return ['None']
    elif context_num == 1:
        return ['A', 'C', 'T', 'G']
    elif context_num == 1.5:
        return ['C*pG', 'CpG*', 'TpC*', 'G*pA',
                'A', 'C', 'T', 'G']
    elif context_num == 2:
        dinucs = [''.join(dinuc)
                  for dinuc in it.combinations_with_replacement('ACTG', 2)]
        return dinucs
    elif context_num == 3:
        trinucs = [''.join(trinuc)
                   for trinuc in it.combinations_with_replacement('ACTG', 3)]
        return trinucs


def compute_mutation_context(bed, gs, df, opts):
    # prepare info for running permutation test
    gene_mut = df[df['Gene']==bed.gene_name]
    cols = ['Chromosome', 'Start_Position', 'Reference_Allele',
            'Tumor_Allele', 'Variant_Classification', 'Protein_Change']
    mut_info = gene_mut[cols]
    gs.set_gene(bed)

    # get sequence context
    if 'seed' in opts:
        sc = permutation2020.python.sequence_context.SequenceContext(gs, seed=opts['seed'])
    else:
        sc = permutation2020.python.sequence_context.SequenceContext(gs)

    # count total mutations in gene
    total_mut = len(mut_info)

    # fix nucleotide letter if gene is on - strand
    if bed.strand == '-':
        mut_info['Tumor_Allele'] = mut_info['Tumor_Allele'].map(lambda x: utils.rev_comp(x))

    # get coding positions, mutations unmapped to the reference tx will have
    # NA for a coding position
    pos_list = []
    for ix, row in mut_info.iterrows():
        coding_pos = bed.query_position(bed.strand, row['Chromosome'], row['Start_Position'])
        pos_list.append(coding_pos)
    mut_info['Coding Position'] = pos_list

    # recover mutations that could not be mapped to the reference transcript
    # for a gene before being dropped (next step)
    unmapped_mut_info = recover_unmapped_mut_info(mut_info, bed, sc, opts)

    # drop mutations wich do not map to reference tx
    mut_info = mut_info.dropna(subset=['Coding Position'])  # mutations need to map to tx
    mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
    unmapped_muts = total_mut - len(mut_info)

    cols = ['Context', 'Tumor_Allele', 'Coding Position']
    if len(mut_info) > 0:
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        mut_info['Context'] = mut_info['Coding Position'].apply(lambda x: sc.pos2context[x])

        # group mutations by context
        unmapped_mut_df = pd.DataFrame(unmapped_mut_info)
        rename_dict = {'Codon Pos': 'Coding Position'}
        unmapped_mut_df = unmapped_mut_df.rename(columns=rename_dict)
        tmp_df = pd.concat([mut_info[cols], unmapped_mut_df[cols]])
        context_cts = tmp_df['Context'].value_counts()
        context_to_mutations = dict((name, group['Tumor_Allele'])
                                    for name, group in tmp_df.groupby('Context'))
    else:
        # initialize empty results if there are no mutations
        context_cts = pd.Series([])
        context_to_mutations = {}
        tmp_df = pd.DataFrame(columns=cols)

    return context_cts, context_to_mutations, tmp_df, gs, sc


def get_chasm_context(tri_nuc):
    """Returns the mutation context acording to CHASM.

    For more information about CHASM's mutation context, look
    at http://wiki.chasmsoftware.org/index.php/CHASM_Overview.
    Essentially CHASM uses a few specified di-nucleotide contexts
    followed by single nucleotide context.

    Parameters
    ----------
    tri_nuc : str
        three nucleotide string with mutated base in the middle.

    Returns
    -------
    chasm context : str
        a string representing the context used in CHASM
    """
    # check if string is correct length
    if len(tri_nuc) != 3:
        raise ValueError('Chasm context requires a three nucleotide string '
                         '(Provided: "{0}")'.format(tri_nuc))

    # try dinuc context if found
    if tri_nuc[1:] == 'CG':
        return 'C*pG'
    elif tri_nuc[:2] == 'CG':
        return 'CpG*'
    elif tri_nuc[:2] == 'TC':
        return 'TpC*'
    elif tri_nuc[1:] == 'GA':
        return 'G*pA'
    else:
        # just return single nuc context
        return tri_nuc[1]


def get_context(chr, pos_list, strand, fa, context_type):
    nuc_contexts = []
    if context_type in [1, 2]:
        # case where context matters
        index_context = int(context_type) - 1  # subtract 1 since python is zero-based index
        for pos in pos_list:
            nucs = fa.fetch(reference=chr,
                            start=pos-index_context,
                            end=pos+1).upper()
            #commented out always using positive strand
            #to use positive strand context uncomment the following
            #if strand == '-':
                #nucs = utils.rev_comp(nucs)
            if 'N' not in nucs:
                nuc_contexts.append(nucs)
            else:
                nuc_contexts.append(None)
    elif context_type in [1.5, 3]:
        # use the nucleotide context from chasm if nuc
        # context is 1.5 otherwise always use a three
        # nucleotide context
        for pos in pos_list:
            nucs = fa.fetch(reference=chr,
                            start=pos-1,
                            end=pos+2).upper()
            #commented out always using positive strand
            #to use positive strand context uncomment the following
            #if strand == '-':
                #nucs = utils.rev_comp(nucs)
            if context_type == 1.5 and nucs:
                nucs = get_chasm_context(nucs)

            if 'N' not in nucs:
                nuc_contexts.append(nucs)
            else:
                nuc_contexts.append(None)
    else:
        nuc_contexts = ['None'] * len(pos_list)

    return nuc_contexts


def get_aa_mut_info(coding_pos, somatic_base, gene_seq):
    """Retrieves relevant information about the effect of a somatic
    SNV on the amino acid of a gene.

    Information includes the germline codon, somatic codon, codon
    position, germline AA, and somatic AA.

    Parameters
    ----------
    coding_pos : iterable of ints
        Contains the base position (0-based) of the mutations
    somatic_base : list of str
        Contains the somatic nucleotide for the mutations
    gene_seq : GeneSequence
        gene sequence

    Returns
    -------
    aa_info : dict
        information about the somatic mutation effect on AA's
    """
    # get codon information into three lists
    ref_codon, codon_pos, pos_in_codon, ref_nuc = it.izip(*[cutils.pos_to_codon(gene_seq, p)
                                                            for p in coding_pos])
    ref_codon, codon_pos, pos_in_codon, ref_nuc = list(ref_codon), list(codon_pos), list(pos_in_codon), list(ref_nuc)

    # construct codons for mutations
    mut_codon = [(list(x) if x != 'Splice_Site' else []) for x in ref_codon]
    for i in range(len(mut_codon)):
        # splice site mutations are not in a codon, so skip such mutations to
        # prevent an error
        if pos_in_codon[i] is not None:
            pc = pos_in_codon[i]
            mut_codon[i][pc] = somatic_base[i]
    mut_codon = [(''.join(x) if x else 'Splice_Site') for x in mut_codon]

    # output resulting info
    aa_info = {'Reference Codon': ref_codon,
               'Somatic Codon': mut_codon,
               'Codon Pos': codon_pos,
               'Reference Nuc': ref_nuc,
               'Reference AA': [(utils.codon_table[r] if (r in utils.codon_table) else None)
                                for r in ref_codon],
               'Somatic AA': [(utils.codon_table[s] if (s in utils.codon_table) else None)
                              for s in mut_codon]}

    return aa_info


def get_unmapped_aa_mut_info(mut_info, genome_fa, strand, chr, context_type):

    # get information on the nucleotide context
    mycontexts = get_context(chr, mut_info['Start_Position'],
                             strand, genome_fa, context_type)

    # get information about the effect of the protein change
    if len(mut_info) > 0:
        not_splice_site = mut_info['Variant_Classification'].map(lambda x: x!='Splice_Site')
        prot_change = [AminoAcid(p) for p in mut_info[not_splice_site]['Protein_Change']]
    else:
        prot_change = []
    codon_pos, germ_aa, somatic_aa = [], [], []
    LARGE_NUMBER = 100000  # sufficiently large number to prevent accidental overlap of codon positions
    tmp_index = 0
    bad_mut_ix, good_mut_ix = [], []
    for i in range(len(mut_info)):
        if not mycontexts[i]:
            bad_mut_ix.append(i)  # remove invalid/missing mutation
            codon_pos.append(None)
            germ_aa.append(None)
            somatic_aa.append(None)
            if not_splice_site.iloc[i]:
                tmp_index += 1
        elif not_splice_site.iloc[i]:
            if prot_change and (not prot_change[tmp_index].is_valid or \
                                prot_change[tmp_index].is_missing_info):
                bad_mut_ix.append(i)  # remove invalid/missing mutation
                codon_pos.append(None)
                germ_aa.append(None)
                somatic_aa.append(None)
                tmp_index += 1
            else:
                good_mut_ix.append(i)
                codon_pos.append(LARGE_NUMBER + prot_change[tmp_index].pos)
                germ_aa.append(prot_change[tmp_index].initial)
                somatic_aa.append(prot_change[tmp_index].mutated)
                tmp_index += 1
        else:
            good_mut_ix.append(i)
            codon_pos.append('Splice_Site')
            germ_aa.append('Splice_Site')
            somatic_aa.append('Splice_Site')

    # remove bad mutations from results
    mycontexts = utils.filter_list(mycontexts, bad_mut_ix)
    germ_aa = utils.filter_list(germ_aa, bad_mut_ix)
    somatic_aa = utils.filter_list(somatic_aa, bad_mut_ix)
    codon_pos = utils.filter_list(codon_pos, bad_mut_ix)

    # information about the effect of mutations that could not be mapped
    # to the reference isoform of a gene.
    tumor_allele = mut_info.iloc[np.array(good_mut_ix, dtype=int)]['Tumor_Allele'].tolist()
    aa_info = {'Context': mycontexts,
               'Codon Pos': codon_pos,
               'Reference AA': germ_aa,
               'Somatic AA': somatic_aa,
               'Tumor_Allele': tumor_allele}

    return aa_info


def recover_unmapped_mut_info(mut_info, bed, sc, opts):
    # retreive info based on annotated protein effects and genomic coordinates
    has_unmapped_opts = ('use_unmapped' in opts) and ('genome' in opts)
    use_unmapped = opts['use_unmapped'] and opts['genome']
    if has_unmapped and use_unmapped:
        genome_fa = pysam.Fastafile(opts['genome'])
        # try to still use mutations that are not on the reference transcript
        tmp_mut_info = mut_info[mut_info['Coding Position'].isnull()]
        unmapped_mut_info = get_unmapped_aa_mut_info(tmp_mut_info,
                                                     genome_fa,
                                                     bed.strand,
                                                     bed.chrom,
                                                     opts['context'])
        genome_fa.close()

        # filter out cases where the nucleotide context does not exist
        # on the reference transcript
        bad_contexts = [i for i in range(len(unmapped_mut_info['Context']))
                        if not sc.is_valid_context(unmapped_mut_info['Context'][i])]
        for key in unmapped_mut_info:
            unmapped_mut_info[key] = utils.filter_list(unmapped_mut_info[key],
                                                       bad_contexts)
    else:
        unmapped_mut_info = {'Context': [], 'Reference AA': [], 'Codon Pos': [],
                             'Somatic AA': [], 'Tumor_Allele': []}
    return unmapped_mut_info
