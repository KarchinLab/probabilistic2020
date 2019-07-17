import numpy as np
import prob2020.python.utils as utils
import prob2020.python.mutation_context


class SequenceContext(object):
    """The SequenceContext class allows for deciphering sequence context
    and for randomly permuting mutation positions while respecting sequence context.
    """

    def __init__(self, gene_seq, seed=None):
        self._init_context(gene_seq)
        self.seed = seed  # seed for random number generator
        context_names = prob2020.python.mutation_context.get_all_context_names(gene_seq.nuc_context)
        self.prng_dict = {
            c: np.random.RandomState(seed=self.seed)
            for c in context_names
        }
        self.prng_dict['N'] = np.random.RandomState(seed=self.seed)

    def _init_context(self, gene_seq):
        """Initializes attributes defining mutation contexts and their position.

        The self.context2pos and self.pos2context dictionaries map from
        sequence context to sequence position and sequence position to
        sequence context, respectively. These attributes allow for randomly
        sampling of mutation positions while respecting sequence context in the
        randomization-based test.

        Parameters
        ----------
        gene_seq : GeneSequence
            GeneSequence object from the gene_sequence module
        """
        self.context2pos, self.pos2context = {}, {}
        gene_len = len(gene_seq.exon_seq)  # get length of CDS
        five_ss_len = 2*len(gene_seq.five_prime_seq)  # total length of 5' splice sites
        three_ss_len = 2*len(gene_seq.three_prime_seq)  # total length of 3' splice sites

        if gene_seq.nuc_context in [1, 2]:
            # case where context matters
            index_context = int(gene_seq.nuc_context) - 1  # subtract 1 since python is zero-based index
            for i in range(index_context, gene_len):
                nucs = gene_seq.exon_seq[i-index_context:i+1]
                self.context2pos.setdefault(nucs, [])
                self.context2pos[nucs].append(i)
                self.pos2context[i] = nucs

            # sequence context for five prime splice site
            for i, five_ss in enumerate(gene_seq.five_prime_seq):
                first_nucs = five_ss[1-index_context:1+1]
                second_nucs = five_ss[2-index_context:2+1]
                first_pos = 2*i + gene_len
                second_pos = 2*i + gene_len + 1
                self.context2pos.setdefault(first_nucs, [])
                self.context2pos[first_nucs].append(first_pos)
                self.context2pos.setdefault(second_nucs, [])
                self.context2pos[second_nucs].append(second_pos)
                self.pos2context[first_pos] = first_nucs
                self.pos2context[second_pos] = second_nucs
            # sequence context for three prime splice site
            for i, three_ss in enumerate(gene_seq.three_prime_seq):
                first_nucs = three_ss[1-index_context:1+1]
                second_nucs = three_ss[2-index_context:2+1]
                first_pos = 2*i + gene_len + five_ss_len
                second_pos = 2*i + gene_len + five_ss_len + 1
                self.context2pos.setdefault(first_nucs, [])
                self.context2pos[first_nucs].append(first_pos)
                self.context2pos.setdefault(second_nucs, [])
                self.context2pos[second_nucs].append(second_pos)
                self.pos2context[first_pos] = first_nucs
                self.pos2context[second_pos] = second_nucs

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
                    context = prob2020.python.mutation_context.get_chasm_context(nucs)
                else:
                    context = nucs
                self.context2pos.setdefault(context, [])
                self.context2pos[context].append(i)
                self.pos2context[i] = context

            # sequence context for five prime splice site
            for i, five_ss in enumerate(gene_seq.five_prime_seq):
                first_nucs = five_ss[:3]
                second_nucs = five_ss[1:4]
                first_pos = 2*i + gene_len
                second_pos = 2*i + gene_len + 1
                if ncontext == 1.5:
                    first_context = prob2020.python.mutation_context.get_chasm_context(first_nucs)
                    second_context = prob2020.python.mutation_context.get_chasm_context(second_nucs)
                else:
                    first_context = first_nucs
                    second_context = second_nucs
                self.context2pos.setdefault(first_context, [])
                self.context2pos[first_context].append(first_pos)
                self.context2pos.setdefault(second_context, [])
                self.context2pos[second_context].append(second_pos)
                self.pos2context[first_pos] = first_context
                self.pos2context[second_pos] = second_context
            # sequence context for three prime splice site
            for i, three_ss in enumerate(gene_seq.three_prime_seq):
                first_nucs = three_ss[:3]
                second_nucs = three_ss[1:4]
                first_pos = 2*i + gene_len + five_ss_len
                second_pos = 2*i + gene_len + five_ss_len + 1
                if ncontext == 1.5:
                    first_context = prob2020.python.mutation_context.get_chasm_context(first_nucs)
                    second_context = prob2020.python.mutation_context.get_chasm_context(second_nucs)
                else:
                    first_context = first_nucs
                    second_context = second_nucs
                self.context2pos.setdefault(first_context, [])
                self.context2pos[first_context].append(first_pos)
                self.context2pos.setdefault(second_context, [])
                self.context2pos[second_context].append(second_pos)
                self.pos2context[first_pos] = first_context
                self.pos2context[second_pos] = second_context

            # hack solution for context for first nuc
            if gene_seq.exon_seq:
                first_nuc = gene_seq.exon_seq[0] + gene_seq.exon_seq[:2]
                if ncontext == 1.5:
                    first_context = prob2020.python.mutation_context.get_chasm_context(first_nuc)
                else:
                    first_context = first_nuc
                self.pos2context[0] = first_context
                self.context2pos.setdefault(first_context, [])
                self.context2pos[first_context].append(0)
                last_nuc = gene_seq.exon_seq[-2:] + gene_seq.exon_seq[-1]
                if ncontext == 1.5:
                    last_context = prob2020.python.mutation_context.get_chasm_context(last_nuc)
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
            for i in range(gene_len + five_ss_len + three_ss_len):
                self.pos2context[i] = 'None'
            self.context2pos['None'] = range(gene_len + five_ss_len + three_ss_len)

    def is_valid_context(self, ctxt):
        """Checks if provided context is valid (previously seen).

        Parameters
        ----------
        ctxt : str
            mutation context
        """
        return ctxt in self.context2pos

    def random_context_pos(self, num, num_permutations, context):
        """Samples with replacement available positions matching the
        sequence context.

        Note: this method does random sampling only for an individual
        sequence context.

        Parameters
        ----------
        num : int
            Number of positions to sample for each permutation. This
            is the number of actually observed mutations having the
            matching sequence context for this gene.
        num_permutations : int
            Number of permutations for permutation test.
        context : str
            Sequence context.

        Returns
        -------
        random_pos : np.array
            num_permutations X num sized array that represents the
            randomly sampled positions for a specific context.
        """
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
        random_pos = self.prng_dict[context].choice(available_pos, (num_permutations, num))
        return random_pos

    def random_pos(self, context_iterable, num_permutations):
        """Obtains random positions w/ replacement which match sequence context.

        Parameters
        ----------
        context_iterable: iterable containing two element tuple
            Records number of mutations in each context. context_iterable
            should be something like [('AA', 5), ...].
        num_permutations : int
            Number of permutations used in the permutation test.

        Returns
        -------
        position_list : list
            Contains context string and the randomly chosen positions
            for that context.
        """
        position_list = []
        for contxt, n in context_iterable:
            pos_array = self.random_context_pos(n, num_permutations, contxt)
            position_list.append([contxt, pos_array])
        return position_list

