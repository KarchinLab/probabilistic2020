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

Running the statistical test
----------------------------

The statistical tests account for gene sequence and mutational context.
Each gene is represented by a single reference transcript (above is longest CDS SNVBox transcript).
By default the relevant sequence context for mutations are utilized from
chasm paper (denoted by **-c 1.5** parameter). This includes some common dinucletoide contexts
like CpG, and otherwise just a single base.

**Technical detail:** Running on the obtained pan-cancer data may take several hours to run on a single
core. Specifying the **-p** parameter to use multiple processors will speed up run time if available.

Running oncogene sub-command
++++++++++++++++++++++++++++

The oncogene sub-command examines missense position clustering (by codon) and elevated
*in silico* pathogenicity scores (VEST).


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

.. code-block:: bash

   $ probabilistic2020 tsg \
        -i genes.fa \
        -b genes.bed \
        -p 10 \
        -c 1.5 \
        -o tsg_output.txt

