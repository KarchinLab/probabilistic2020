# 20/20 Permutation Test

The 20/20 permutation test identifies signficant oncogenes and tumor suppressor genes 
(TSG). Putative signficant oncogenes are found through evaluating the position of 
missense mutations (clustered missense mutations tend to indicate actiavting mutations).
While statistically signficant TSGs are found by abnormally high number of deleterious
mutations.

The 20/20 Permutation Test performs a permutation test by matching
null mutation context with the observed mutation context. As such,
the permutation test avoids using a background distribution for genes.
Significance is determined by considering only mutations within the
same gene.

The 20/20 Permutation Test has nice properties since it accounts
for several factors that could effect the significance of driver genes.

* gene length
* mutation context
* codon bias

## Creating the gene annotation

To perform permutations, the gene structure (ie position of exon, coding region) needs to 
be defined. For each gene, only one transcript is used. The list of genes utilized and
their corresponding transcript is obtained through [SNVBox](http://wiki.chasmsoftware.org/index.php/Main_Page). Follow instructions on installing SNVBox, if you want to recreate this
entire process. Assuming SNVBox is installed, then the longest transcript is selected
for each gene in SNVBox:

    $ mysql [options] < scripts/longest_snvbox_tx.sql > output/gene_bed/longest_tx.txt

The longest_snvbox_tx.sql script is also located on [github](https://gist.github.com/ctokheim/18363041037e375f411c). 
However if multiple transcripts have the same maximum amino acid length, then multiple 
transcripts will be reported in the longest_tx.txt file. To remove redundancies 
arbitrarily by selecting the first transcript occurrence, the scripts/unique_tx.py script 
is used:

    $ python scripts/unique_tx.py -i output/gene_bed/longest_tx.txt -o output/gene_bed/unique_longest_tx.txt

A list of transcript IDs are obtained using `cut` for input into the UCSC table browser:

    $ cut -f2,2 output/gene_bed/unique_longest_tx.txt | tail -n +2 > output/gene_bed/tx_ids.txt

Then obtain a BED file from the table browser using both the refGene and ensGene tracks.
The output was saved as output/gene_bed/ucsc_table_browser_output.bed. The output from the UCSC table browser is then associated with a gene name rather than a transcript ID by the following:

    $ python scripts/create_gene_bed.py -b output/gene_bed/ucsc_table_browser_output.bed -g output/gene_bed/unique_longest_tx.txt -o data/snvboxGenes.bed

NOTE: Any transcripts from the table browser that had multiple positions (eg on chrX 
and chrY) along with their corresponding gene were not saved in the final output.

## Creating Gene FASTA

    $ python bin/extract_gene_seq.py --log-level=DEBUG -i data/hg19.fa -b data/snvboxGenes.bed -o data/snvboxGenes.fa
