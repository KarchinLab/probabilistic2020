Advanced 
========

If you need to run on different genomes (e.g. mouse) or different reference transcripts, this will
require creating correct input files.

* gene annotations (BED file) 
* gene sequence (FASTA file)
* conservation/vest scores
* Needed for 20/20+

  * generating gene features (observed and simulated)
  * generating simulated MAF files


Creating the gene annotation (BED file) from SNVBox
---------------------------------------------------

To perform permutations, the gene structure (ie position of exon, coding region) needs to 
be defined. For each gene, only one transcript is used. The list of genes utilized and
their corresponding transcript is obtained through `SNVBox <http://wiki.chasmsoftware.org/index.php/Main_Page>`_. Follow instructions on installing SNVBox, if you want to recreate this
entire process. Assuming SNVBox is installed, then the longest transcript is selected
for each gene in SNVBox:

.. code-block:: bash

   $ mysql [options] < scripts/longest_snvbox_tx.sql > output/gene_bed/longest_tx.txt

The longest_snvbox_tx.sql script is also located on `github <https://gist.github.com/ctokheim/18363041037e375f411c>`_. 
However if multiple transcripts have the same maximum amino acid length, then multiple 
transcripts will be reported in the longest_tx.txt file. To remove redundancies 
arbitrarily by selecting the first transcript occurrence, the scripts/unique_tx.py script 
is used:

.. code-block:: bash

   $ python scripts/unique_tx.py -i output/gene_bed/longest_tx.txt -o output/gene_bed/unique_longest_tx.txt

A list of transcript IDs are obtained using `cut` for input into the UCSC table browser:

.. code-block:: bash

   $ cut -f2,2 output/gene_bed/unique_longest_tx.txt | tail -n +2 > output/gene_bed/tx_ids.txt

Then obtain a BED file from the table browser using both the refGene and ensGene tracks.
The output was saved as output/gene_bed/ucsc_table_browser_output.bed. The output from the UCSC table browser is then associated with a gene name rather than a transcript ID by the following:

.. code-block:: bash

   $ python scripts/create_gene_bed.py -b output/gene_bed/ucsc_table_browser_output.bed -i output/gene_bed/ignore_chroms.txt -g output/gene_bed/unique_longest_tx.txt -o data/snvboxGenes.bed

NOTE: Any transcripts from the table browser that had multiple positions (eg on chrX 
and chrY) along with their corresponding gene were not saved in the final output. To prevent genes from not being used because of duplicate positions but on say non-canonical chromosome names, use the -i option to filter out those issues.

Creating your own Gene FASTA
----------------------------

Gene sequences are extracted from a genome FASTA file. To do this, you need
a BED file with names corresponding to genes, and a genome FASTA (e.g. hg19).
Obtaining a BED file can be done using the above steps. Creating the gene
sequence FASTA is then done by the `extract_gene_seq` script:

.. code-block:: bash

    $ extract_gene_seq -i mygenome.fa -b mygenes.bed -o mygenes.fa

Where mygenome.fa is the genome FASTA file, mygenes.bed contains a single reference transcript for each gene (with gene names not transcript names), and mygenes.fa is the FASTA
file separated out by genes.

Extracting conservation and VEST scores
---------------------------------------

Needed for 20/20+
-----------------

Generating gene features
++++++++++++++++++++++++

Generating simulated MAF file
+++++++++++++++++++++++++++++
