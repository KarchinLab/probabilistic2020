"""Fetches gene sequence from gene fasta created by extract_genes.py"""


class GeneSequence(object):

    def __init__(self, fasta_obj,
                 nuc_context=1):
        self.fasta = fasta_obj
        self.nuc_context = nuc_context

    def set_gene(self, bed_line):
        self.bed = bed_line  # gene that was specified as BED
        self._reset_seq()  # fetch sequence for bed line

    def _reset_seq(self):
        exon_seq_list, five_ss_seq_list, three_ss_seq_list = self._fetch_seq()
        self.exon_seq = ''.join(exon_seq_list)
        self.three_prime_seq = [(exon_seq_list[i][-1:] + ss if self.nuc_context == 2 else ss)
                                for i, ss in enumerate(three_ss_seq_list)]
        self.five_prime_seq = [(exon_seq_list[i+1][:1] + ss if self.nuc_context == 2 else ss)
                               for i, ss in enumerate(three_ss_seq_list)]
        self._to_upper()  # make sure all sequences are in upper case

    def _to_upper(self):
        self.exon_seq = self.exon_seq.upper()
        self.three_prime_seq = [s.upper() for s in self.three_prime_seq]
        self.five_prime_seq = [s.upper() for s in self.five_prime_seq]

    def _fetch_seq(self):
        exons = []
        three_prime_ss = []
        five_prime_ss = []
        num_exons = self.bed.get_num_exons()
        for i in range(num_exons):
            # add exon sequence
            tmp_id = '{0};exon{1}'.format(self.bed.gene_name, i)
            tmp_exon = self.fasta.fetch(reference=tmp_id)
            exons.append(tmp_exon)

            # add splice site sequence
            tmp_id_3ss = '{0};3SS'.format(tmp_id)
            tmp_id_5ss = '{0};3SS'.format(tmp_id)
            if num_exons == 1:
                pass
            elif i == 0:
                tmp_3ss = self.fasta.fetch(tmp_id_3ss)
                three_prime_ss.append(tmp_3ss)
            elif i == (num_exons - 1):
                tmp_5ss = self.fasta.fetch(tmp_id_5ss)
                five_prime_ss.append(tmp_5ss)
            else:
                tmp_3ss = self.fasta.fetch(tmp_id_3ss)
                tmp_5ss = self.fasta.fetch(tmp_id_5ss)
                three_prime_ss.append(tmp_3ss)
                five_prime_ss.append(tmp_5ss)
        return exons, five_prime_ss, three_prime_ss

    def pos_to_codon(self, pos):
        codon_pos = pos // 3
        codon_start = codon_pos * 3
        codon = self.exon_seq[codon_start:codon_start+3]
        pos_in_codon = pos % 3
        return codon, codon_pos, pos_in_codon
