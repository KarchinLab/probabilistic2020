"""Fetches gene sequence from gene fasta created by extract_genes.py"""
import prob2020.python.utils as utils


class GeneSequence(object):

    def __init__(self, fasta_obj,
                 nuc_context=1.5):
        self.fasta = fasta_obj
        self.nuc_context = nuc_context

    def set_gene(self, bed_line):
        """Updates gene sequence for a new gene (bed line).

        Parameters
        ----------
        bed_line : BedLine
            BedLine object representing a single gene in a BED file
        """
        self.bed = bed_line  # gene that was specified as BED
        self._reset_seq()  # fetch sequence for bed line

    def _reset_seq(self):
        """Updates attributes for gene represented in the self.bed attribute.

        Sequences are always upper case.
        """
        exon_seq_list, five_ss_seq_list, three_ss_seq_list = self._fetch_seq()
        self.exon_seq = ''.join(exon_seq_list)
        self.three_prime_seq = three_ss_seq_list
        self.five_prime_seq =  five_ss_seq_list
        self._to_upper()  # make sure all sequences are in upper case

    def add_germline_variants(self, germline_nucs, coding_pos):
        """Add potential germline variants into the nucleotide sequence.

        Sequenced individuals may potentially have a SNP at a somatic mutation position.
        Therefore they may differ from the reference genome. This method updates the gene
        germline gene sequence to match the actual individual.

        Parameters
        ----------
        germline_nucs : list of str
            list of DNA nucleotides containing the germline letter
        coding_pos : int
            0-based nucleotide position in coding sequence

        NOTE: the self.exon_seq attribute is updated, no return value
        """
        if len(germline_nucs) != len(coding_pos):
            raise ValueError('Each germline nucleotide should have a coding position')

        es = list(self.exon_seq)
        for i in range(len(germline_nucs)):
            gl_nuc, cpos = germline_nucs[i].upper(), coding_pos[i]
            if not utils.is_valid_nuc(gl_nuc):
                raise ValueError('{0} is not a valid nucleotide'.format(gl_nuc))
            if cpos >= 0:
                es[cpos] = gl_nuc
        self.exon_seq = ''.join(es)

    def _to_upper(self):
        """Convert sequences to upper case."""
        self.exon_seq = self.exon_seq.upper()
        self.three_prime_seq = [s.upper() for s in self.three_prime_seq]
        self.five_prime_seq = [s.upper() for s in self.five_prime_seq]

    def _fetch_seq(self):
        """Fetches gene sequence from PySAM fasta object.

        Returns
        -------
        exons : list of str
            list of exon nucleotide sequences
        five_prime_ss : list of str
            list of 5' splice site sequences
        three_prime_ss : list of str
            list of 3' splice site sequences
        """
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
            tmp_id_5ss = '{0};5SS'.format(tmp_id)
            if num_exons == 1:
                pass
            elif i == 0:
                tmp_5ss = self.fasta.fetch(tmp_id_5ss)
                five_prime_ss.append(tmp_5ss)
            elif i == (num_exons - 1):
                tmp_3ss = self.fasta.fetch(tmp_id_3ss)
                three_prime_ss.append(tmp_3ss)
            else:
                tmp_3ss = self.fasta.fetch(tmp_id_3ss)
                tmp_5ss = self.fasta.fetch(tmp_id_5ss)
                three_prime_ss.append(tmp_3ss)
                five_prime_ss.append(tmp_5ss)
        return exons, five_prime_ss, three_prime_ss


def _fetch_5ss_fasta(fasta, gene_name, exon_num,
                     chrom, strand, start, end):
    """Retreives the 5' SS sequence flanking the specified exon.

    Returns a string in fasta format with the first line containing
    a ">" and the second line contains the two base pairs of 5' SS.

    Parameters
    ----------
    fasta : pysam.Fastafile
        fasta object from pysam
    gene_name : str
        gene name used for fasta seq id
    exon_num : int
        the `exon_num` exon, used for seq id
    chrom : str
        chromsome
    strand : str
        strand, {'+', '-'}
    start : int
        0-based start position
    end : int
        0-based end position

    Returns
    -------
    ss_fasta : str
        string in fasta format with first line being seq id
    """
    if strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end-1,
                             end=end+3)
    elif strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-3,
                             end=start+1)
        ss_seq = utils.rev_comp(ss_seq)
    ss_fasta = '>{0};exon{1};5SS\n{2}\n'.format(gene_name,
                                                exon_num,
                                                ss_seq.upper())
    return ss_fasta


def _fetch_3ss_fasta(fasta, gene_name, exon_num,
                     chrom, strand, start, end):
    """Retreives the 3' SS sequence flanking the specified exon.

    Returns a string in fasta format with the first line containing
    a ">" and the second line contains the two base pairs of 3' SS.

    Parameters
    ----------
    fasta : pysam.Fastafile
        fasta object from pysam
    gene_name : str
        gene name used for fasta seq id
    exon_num : int
        the `exon_num` exon, used for seq id
    chrom : str
        chromsome
    strand : str
        strand, {'+', '-'}
    start : int
        0-based start position
    end : int
        0-based end position

    Returns
    -------
    ss_fasta : str
        string in fasta format with first line being seq id
    """

    if strand == '-':
        ss_seq = fasta.fetch(reference=chrom,
                             start=end-1,
                             end=end+3)
        ss_seq = utils.rev_comp(ss_seq)
    elif strand == '+':
        ss_seq = fasta.fetch(reference=chrom,
                             start=start-3,
                             end=start+1)
    ss_fasta = '>{0};exon{1};3SS\n{2}\n'.format(gene_name,
                                                exon_num,
                                                ss_seq.upper())
    return ss_fasta


def fetch_gene_fasta(gene_bed, fasta_obj):
    """Retreive gene sequences in FASTA format.

    Parameters
    ----------
    gene_bed : BedLine
        BedLine object representing a single gene
    fasta_obj : pysam.Fastafile
        fasta object for index retreival of sequence

    Returns
    -------
    gene_fasta : str
        sequence of gene in FASTA format
    """
    gene_fasta = ''
    strand = gene_bed.strand
    exons = gene_bed.get_exons()
    if strand == '-':
        exons.reverse()  # order exons 5' to 3', so reverse if '-' strand

    # iterate over exons
    for i, exon in enumerate(exons):
        exon_seq = fasta_obj.fetch(reference=gene_bed.chrom,
                                   start=exon[0],
                                   end=exon[1]).upper()
        if strand == '-':
            exon_seq = utils.rev_comp(exon_seq)
        exon_fasta = '>{0};exon{1}\n{2}\n'.format(gene_bed.gene_name,
                                                  i, exon_seq)

        # get splice site sequence
        if len(exons) == 1:
            # splice sites don't matter if there is no splicing
            ss_fasta = ''
        elif i == 0:
            # first exon only, get 3' SS
            ss_fasta = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1])
        elif i == (len(exons) - 1):
            # last exon only, get 5' SS
            ss_fasta = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                        gene_bed.chrom, strand, exon[0], exon[1])
        else:
            # middle exon, get bot 5' and 3' SS
            fasta_3ss = _fetch_3ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1])
            fasta_5ss = _fetch_5ss_fasta(fasta_obj, gene_bed.gene_name, i,
                                         gene_bed.chrom, strand, exon[0], exon[1])
            ss_fasta = fasta_5ss + fasta_3ss

        gene_fasta += exon_fasta + ss_fasta

    return gene_fasta
