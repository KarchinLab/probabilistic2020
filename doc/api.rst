API
===

If you follow the installation instructions for installing as a python package
then the simulation code can be repurposed for other tasks.

Mutation Simulations
--------------------

The first step is to collect a set of annotations include a BED file, gene sequence
FASTA (obtained from extract_gene_seq.py), and a set of mutations in the standard
format for probabilistic 20/20 (cite annotation page). Once those are prepared,
they may be used as input for creating simulated mutations. For this exmaple
I will assume the files are named as follows:

* mutations.txt
* genes.bed
* genes.fa

Gene annotation (BED file)
++++++++++++++++++++++++++

BED files are read into a dictionary with keys as chromosome names and values
as lists of python BED objects.

.. code-block:: python

   import prob2020.python.utils as utils

   # read in BED file
   bed_dict = utils.read_bed('genes.bed')

To get a single BED object, you must access a single element in a list.

.. code-block:: python

   bed_obj = bed_dict['chr1'][0]  # first read gene in chromosome 1

Although in this tutorial a specific gene annotation is used, you may
rather want to iterate over each bed annotation found in bed_dict.

Gene Sequence (FASTA file)
++++++++++++++++++++++++++

Reading in the gene sequence requires `pysam <http://pysam.readthedocs.org/en/latest/api.html>`_ which should already be installed after installation.
The mutational context will also be specified when creating a gene sequence
object. This context follows the same convention as the command line utilities.

.. code-block:: python

   # required imports
   import pysam
   from prob2020.python.gene_sequence import GeneSequence

   # get gene sequence object
   gene_fa = pysam.Fastafile('genes.fa')
   gene_seq = GeneSequence(gene_fa, nuc_context=1.5)

Mutations file
++++++++++++++

The mutations are read by using `pandas <http://pandas.pydata.org/>`_ which is a 
required dependency of this package.

.. code-block:: python

   # required import
   import pandas as pd

   # read in mutations
   mut_df = pd.read_csv('mutations.txt', sep='\t')

Performing simulation
+++++++++++++++++++++

.. code-block:: python

   # required imports
   import prob2020.python.mutation_context as mc
   import prob2020.python.permutation as pm

   # specify options
   opts = {
       'seed': None,  # use random seed for simulation
   }

   # determine number of random permutations to perform
   num_permutations = 10  # you may want this 10,000+ for real data

   # process 
   tmp = mc.compute_mutation_context(bed_obj, gene_seq, mut_df, opts)
   context_cts, context_to_mutations, mutations_df, gs, sc = tmp

   # perform simulations
   sim_result = pm.maf_permutation(context_cts, context_to_mutations,
                                   sc, gs, num_permutations)

sim_result is a list of lists. Each list is a single mutation. 
Each mutation is simulated num_permutations number of times.
