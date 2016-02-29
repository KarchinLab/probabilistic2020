Probabilistic 20/20
===================

The Probabibilistic 20/20 test identifies genes with signficant oncogene-like and tumor suppressor gene-like mutational patterns. 
Putative signficant oncogenes are found through evaluating the position of 
missense mutations (clustered missense mutations tend to indicate actiavting mutations).
While statistically signficant tumor suppressor genes (TSGs) are found by abnormally high number of inactivating mutations.

Probabilistic 20/20 evaluates statistical significance by employing 
monte carlo simulations, which incorporates observed mutation context. Monte carlo
simulations are performed within the same gene and thus avoid building a background
distribution based on other genes.  

The Probabilistic 20/20 test has nice properties since it accounts
for several factors that could effect the significance of driver genes.

* gene length
* mutation context
* codon bias

Installation
------------

Downloading Probabilistic 20/20
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source files can be cloned from github using `git <http://git-scm.com/>`_:

.. code-block:: bash

    $ git clone https://github.com/ctokheim/probabilistic2020.git

The source files can also be manually downloaded from github at https://github.com/ctokheim/probabilistic2020.

Probabilstic 20/20 can be installed either locally or into your python distribution as a package. 

Python Package Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the python package installation, all the required python packages for the 20/20 permutation test will automatically be installed for you.

If you are using a system wide python for installation, you will use the following command.

.. code-block:: bash

    $ sudo pip install probabilistic2020-0.1.0.tar.gz

If your python is locally installed then you do not need to use `sudo`.

.. code-block:: bash

    $ pip install probabilistic2020-0.1.0.tar.gz

The scripts for Probabilstic 20/20 can then be found in `Your_Python_Root_Dir/bin`.

Local installation
~~~~~~~~~~~~~~~~~~

Local installation is a good option if you do not have privilege to install a python package and already have the required packages.

**Required packages:**

* numpy
* scipy
* pandas>=0.17.0
* pysam

If you don't have the above required packages, you will need to install them. For the following commands to work you will need `pip <http://pip.readthedocs.org/en/latest/installing.html>`_. If you are using a system wide python, you will need to use `sudo` before the pip command.

.. code-block:: bash

    $ cd probabilstic2020
    $ pip install -r requirements.txt

Next you will need to build the Probabilistic 20/20 source files. This is can be accomplished in one command.

.. code-block:: bash

    $ make build

Once finished building you can then use the scripts in the `probabilstic2020/bin` directory.

Identifying Non-coding Indels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Retreive gene annotations where each line in a list is an exon in BED format.

Next, combine merge gene annotations with simple repeats.

$ cat data/knownGene.bed data/ensembl.bed data/low_complexity_repeats.bed | sort -k1,1 -k2,2n | ~/software/bedtools/bin/mergeBed -i stdin > data/non_coding_black_list.bed

BEDTOOLs mergeBed seems to provide equivalent merging within a single file as bedops:

$ /projects/clonal-evolution/Mouse/src/bedops_suite/bedops --merge data/non_coding_black_list.bed > data/non_coding_black_list.merged.bed  # same "wc -l" length 

Next, gzip the black list file so that it can be indexed by Tabix in pysam

$ gzip data/non_coding_black_list.bed

Then filter out INDELs which occur in the black list

$ python scripts/non_coding_indel.py -i data/lawrence_indels.txt -b data/non_coding_black_list.bed.gz -o data/non_coding_indels.txt

Calculate non-coding indel background rate:

$ python scripts/calc_non_coding_frameshift_rate.py -b data/non_coding_black_list.merged.bed -g ~/software/bedtools/genomes/human.hg19.genome -i data/non_coding_indels.txt -t 10 -bins 10 -o data/non_coding_fs.background.txt 

