import numpy as np


class SequenceContext(object):

    def __init__(self, gene_seq):
        self._init_context(gene_seq)

    def _init_context(self, gene_seq):
        self.context2pos, self.pos2context = {}, {}

        if gene_seq.nuc_context in [1, 2]:
            # case where context matters
            index_context = gene_seq.nuc_context - 1  # subtract 1 since python is zero-based index
            for i in range(index_context, len(gene_seq.exon_seq)):
                nucs = gene_seq.exon_seq[i-index_context:i+1]
                self.context2pos.setdefault(nucs, [])
                self.context2pos[nucs].append(i)
                self.pos2context[i] = nucs

            # hack solution for context for first nuc
            if gene_seq.exon_seq and gene_seq.nuc_context > 1:
                self.pos2context[0] = gene_seq.exon_seq[0] * 2
                self.context2pos.setdefault(gene_seq.exon_seq[0]*2, [])
                self.context2pos[gene_seq.exon_seq[0]*2].append(0)
        elif gene_seq.nuc_context in [1.5, 3]:
            # use the nucleotide context from chasm if nuc
            # context is 1.5 otherwise always use a three
            # nucleotide context
            ncontext = gene_seq.nuc_context
            for i in range(1, len(gene_seq.exon_seq)-1):
                nucs = gene_seq.exon_seq[i-1:i+2]
                if ncontext == 1.5:
                    context = self.get_chasm_context(nucs)
                else:
                    context = nucs
                self.context2pos.setdefault(context, [])
                self.context2pos[context].append(i)
                self.pos2context[i] = context

            # hack solution for context for first nuc
            if gene_seq.exon_seq:
                first_nuc = gene_seq.exon_seq[0] + gene_seq.exon_seq[:2]
                if ncontext == 1.5:
                    first_context = self.get_chasm_context(first_nuc)
                else:
                    first_context = first_nuc
                self.pos2context[0] = first_context
                self.context2pos.setdefault(first_context, [])
                self.context2pos[first_context].append(0)
                last_nuc = gene_seq.exon_seq[-2:] + gene_seq.exon_seq[-1]
                if ncontext == 1.5:
                    last_context = self.get_chasm_context(last_nuc)
                else:
                    last_context = last_nuc
                last_pos = len(gene_seq.exon_seq) - 1
                self.pos2context[last_pos] = first_context
                self.context2pos.setdefault(last_context, [])
                self.context2pos[last_context].append(last_pos)
        else:
            # case where there is no context,
            # mutations occur with uniform probability at each
            # position
            for i in range(len(gene_seq.exon_seq)):
                self.pos2context[i] = 'None'
            self.context2pos['None'] = range(len(gene_seq.exon_seq))

    def get_chasm_context(self, tri_nuc):
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

    def is_valid_context(self, ctxt):
        return ctxt in self.context2pos

    def random_context_pos(self, num, num_permutations, context):
        # make sure provide context is valid
        if not self.is_valid_context(context):
            error_msg = 'Context ({0}) was never seen in sequence.'.format(context)
            raise ValueError(error_msg)

        # make sure sampling is a positive integer
        if num < 1:
            error_msg = ('There must be at least one sample (specified {0}) '
                         'for a context'.format(num))
            raise ValueError(error_msg)

        # randomly select from available positions that fit the specified context
        available_pos = self.context2pos[context]
        prng = np.random.RandomState()
        random_pos = prng.choice(available_pos, (num_permutations, num))
        return random_pos

    def random_pos(self, context_iterable, num_permutations):
        position_list = []
        for contxt, n in context_iterable:
            pos_array = self.random_context_pos(n, num_permutations, contxt)
            position_list.append([contxt, pos_array])
        return position_list

