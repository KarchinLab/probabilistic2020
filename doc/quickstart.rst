Quick Start
===========

This provides a quick start to running probabilistic2020 to setup
the minimum number of steps to execute probabilistic 20/20 statistical test.


Installation
------------

Please see the `installation page <>`_.

Gene BED annotation
-------------------

BED gene annotation files should contain a single reference transcript per gene. 
The name field in the BED file should contain the gene name (not the transcript).
An example BED file containg the annotations for the largest transcripts in SNVBox 
can be obtained `here <>`_. The 

Creating Gene FASTA
-------------------

Gene sequences are extracted from a genome FASTA file. To do this, you need
a BED file with names corresponding to genes, and a genome FASTA (e.g. hg19).
Creating the gene sequence FASTA is then done by the `extract_gene_seq` script:

.. code-block:: bash

    $ extract_gene_seq -i hg19.fa -b snvboxGenes.bed -o snvboxGenes.fa

In this case the BED file is created using SNVBox, a genome FASTA file for hg19 (hg19.fa), and the
resulting coding sequences for the gene are stored in snvboxGenes.fa.

Mutation Annotation Format (MAF) file
-------------------------------------

Mutations are saved in a MAF-like format. Not All fields in MAF spec are required,
and columns may be in any order. You can download a small example of mutations
for colorectal adenocarcinoma `here <>`_. 

Running an Example
------------------

To execute the statistical test for TSG-like genes by examining elevated proportion 
of inactivating mutations, the `tsg` sub-command  for `probabilistic2020` is used.
To limit the run time for this example, you can limit the number of iterations to
10,000 with the `-n` parameter.

.. code-block:: bash

    $ probabilistic2020 tsg \
        -n 10000 \
        -i snvboxgenes.fa \
        -b snvboxGenes.bed \
        -m colorectal.txt \
        -o colorectal_output.txt

