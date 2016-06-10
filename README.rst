Probabilistic 20/20
===================

The Probabibilistic 20/20 test identifies genes with signficant oncogene-like and tumor suppressor gene-like mutational patterns for small coding region variants. 
Putative signficant oncogenes are found through evaluating 
missense mutation clustering and *in silico* pathogenicity scores. Often highly clustered missense
mutations are indicative of activating mutations.
While statistically signficant tumor suppressor genes (TSGs) are found by abnormally high proportion of inactivating mutations.

Probabilistic 20/20 evaluates statistical significance by employing 
monte carlo simulations, which incorporates observed mutation context. Monte carlo
simulations are performed within the same gene and thus avoid building a background
distribution based on other genes. This means that the statistical test can be applied 
to either all genes in the exome from exome sequencing or to a certain target set of genes
from targeted sequencing.

The Probabilistic 20/20 test has nice properties since it accounts
for several factors that could effect the significance of driver genes.

* gene length
* mutation context
* gene sequence (e.g. codon bias)

Documentation
-------------

.. image:: http://readthedocs.org/projects/probabilistic2020/badge/?version=latest
    :target: http://probabilistic2020.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Please see the `documentation <http://probabilistic2020.readthedocs.org/>`_ on readthedocs for more details.


Installation
------------

.. image:: https://travis-ci.org/KarchinLab/probabilistic2020.svg?branch=master
    :target: https://travis-ci.org/KarchinLab/probabilistic2020


Python Package Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the python package installation, all the required python packages for the probabibilistic 20/20 test will automatically be installed for you.

To install the package into python you can use `pip`. If you are installing to a system wide python then you may need to use `sudo` before the pip command.

.. code-block:: bash

    $ pip install probabilistic2020

The scripts for Probabilstic 20/20 can then be found in `Your_Python_Root_Dir/bin`. You can
check the installation with the following:

.. code-block:: bash

    $ which probabilistic2020
    $ probabibilistic2020 --help

Local installation
~~~~~~~~~~~~~~~~~~

Local installation is a good option if you do not have privilege to install a python package and already have the required packages.  The source files can also be manually downloaded from github at https://github.com/KarchinLab/probabilistic2020/releases.

**Required packages:**

* numpy
* scipy
* pandas>=0.17.0
* pysam

If you don't have the above required packages, you will need to install them. For the following commands to work you will need `pip <http://pip.readthedocs.org/en/latest/installing.html>`_. If you are using a system wide python, you will need to use `sudo` before the pip command. Also if you are using python 3.X then you likely will have to install pysam version >=0.9.0.

.. code-block:: bash

    $ cd probabilstic2020
    $ pip install -r requirements.txt

If you want the exact package version used for development on python 2.7, then instead use the requirements_dev.txt. Next you will need to build the Probabilistic 20/20 source files. This is can be accomplished in one command.

.. code-block:: bash

    $ make build

Once finished building you can then use the scripts in the `probabilstic2020/prob2020/console` directory. You can check the build worked by the following:

.. code-block:: bash

    $ python prob2020/cosole/probabibilistic2020.py --help
