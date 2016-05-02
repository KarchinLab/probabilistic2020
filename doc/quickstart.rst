Quick Start
===========

The quick start is meant to test that everything is working with the installation
of the probabilistic2020 package.
This provides running probabilistic2020 with
the minimum number of steps to execute the statistical test.
For more expansive user instructions see :ref:`tutorial-ref`.

Installation
------------

Please see the :ref:`install-ref`.

Downloading Example
-------------------

Download the quick start example data, and extract the resulting tarball.

.. code-block:: bash

    $ wget http://karchinlab.org/data/2020+/pancreatic_example.tar.gz
    $ tar xvzf pancreatic_example.tar.gz
    $ cd pancreatic_example

Input files
-----------

Gene BED annotation
~~~~~~~~~~~~~~~~~~~

BED gene annotation files should contain a single reference transcript per gene. 
The name field in the BED file should contain the gene name (not the transcript).
An example BED file containg the annotations for the largest transcripts in SNVBox 
is named snvboxGenes.bed. 

Gene FASTA
~~~~~~~~~~

Gene sequences are extracted from a genome FASTA file, and is a step that only needs to be done once. This has already been done for the example BED file provided, but if you were to use a different transcript annotation then you would need to follow the :ref:`make-fasta`.

Mutation Annotation Format (MAF) file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mutations are saved in a MAF-like format. Not All fields in MAF spec are required,
and columns may be in any order. Mutations for pancreatic adenocarcinoma are in the
file pancreatic_adenocarcinoma.txt.

Running the Example
-------------------

To execute the statistical test for TSG-like genes by examining elevated proportion 
of inactivating mutations, the **tsg** sub-command  for **probabilistic2020** is used.
To limit the run time for this example, you can limit the number of iterations to
10,000 with the **-n** parameter.

.. code-block:: bash

    $ probabilistic2020 tsg \
        -n 10000 \
        -i snvboxGenes.fa \
        -b snvboxGenes.bed \
        -m pancreatic_adenocarcinoma.txt \
        -o pancreatic_output.txt
