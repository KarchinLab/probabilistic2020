.. _tutorial-ref:

Tutorial
========

Probabilistic 20/20 consists of two broad statistical tests (oncogene-like and tsg-like) 
and somatic mutation simulation framework. Internally, the simulation framework is 
used to establish statistical significance in the hypothesis test through the 
**probabilistic2020** command. However, the simulation framework through the **mut_annotate** command can 
also be used to create a simulated MAF file where aferwords all mutations are distributed
like passengers based on uniform null distribution. Moreover, a set of mutational
features for each gene representative of driver genes (used in 20/20+) can also be
created.

In this tutorial we go through a some what computationally intensive Pan-cancer
analysis. Please see the quick start for a fast example.

Data setup
----------

First, you need to obtain the data necessary to run the pan-cancer
example.

.. code-block:: bash

   $ wget /path/to/fasta
   $ wget /path/to/genes.bed
   $ extract_gene_seq
   $ wget /path/to/pancan.maf
   $ wget /path/to/score_dir
   $ tar xvzf path_to_score_tarball

.. _make-fasta:

Creating Gene FASTA
~~~~~~~~~~~~~~~~~~~

Gene sequences are extracted from a genome FASTA file, and is a step that only needs to be done once.  
To do this, you need a BED file with names corresponding to genes, and a genome FASTA (e.g. hg19).
You can download hg19 from `here <http://karchinlab.org/data/2020+/hg19.fa.gz>`_.
Creating the gene sequence FASTA is then done by the `extract_gene_seq` script:

.. code-block:: bash

    $ extract_gene_seq -i hg19.fa -b snvboxGenes.bed -o snvboxGenes.fa

In this case the BED file is created using SNVBox, a genome FASTA file for hg19 (hg19.fa), and the
resulting coding sequences for the gene are stored in snvboxGenes.fa.

Running the statistical test
----------------------------

The statistical tests account for gene sequence and mutational context.
Each gene is represented by a single reference transcript (above is longest CDS SNVBox transcript).
By default the relevant sequence context for mutations are utilized from
`CHASM paper <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2763410/>`_ (denoted by **-c 1.5** parameter). This includes some common dinucletoide contexts
like CpG, and otherwise just a single base. Ultimately a multiple testing corrected q-value
is reported using the Benjamini-Hochberg (BH) method.

**Technical detail:** Running on the obtained pan-cancer data may take several hours to run on a single
core. Specifying the **-p** parameter to use multiple processors will speed up run time if available.
Lowering the number of iterations (default: 100,000) will decrease run time, but also decrease the resolution
of p-values.

Running oncogene sub-command
++++++++++++++++++++++++++++

The oncogene sub-command examines missense position clustering (by codon) and elevated
*in silico* pathogenicity scores (VEST). The p-values will be combined using fischer's method
to report a single p-value with a BH FDR.

.. code-block:: bash

   $ probabilistic2020 oncogene \
        -i genes.fa \
        -b genes.bed \
        -s score_dir \
        -c 1.5
        -p 10 \
        -o oncogene_output.txt

Running tsg sub-command
+++++++++++++++++++++++

Evaluating for elevated proportion of inactivating point mutations to find TSG-like genes,
can be done using the **tsg** sub-command.

.. code-block:: bash

   $ probabilistic2020 tsg \
        -i genes.fa \
        -b genes.bed \
        -p 10 \
        -c 1.5 \
        -o tsg_output.txt
