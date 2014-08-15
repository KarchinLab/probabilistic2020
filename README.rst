20/20 Permutation Test
======================

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

Installation
------------

Downloading 20/20 Permutation Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source files can be cloned from github using `git <http://git-scm.com/>`_:

.. code-block:: bash

    $ git clone https://github.com/ctokheim/permutation2020.git

The source files can also be manually downloaded from github at https://github.com/ctokheim/permutation2020.

20/20 permutation test can be installed either locally or into your python distribution as a package. 

Python Package Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the python package installation, all the required python packages for the 20/20 permutation test will automatically be installed for you.

If you are using a system wide python for installation, you will use the following command.

.. code-block:: bash

    $ sudo pip install permutation2020-0.1.0.tar.gz

If your python is locally installed then you do not need to use `sudo`.

.. code-block:: bash

    $ pip install permutation2020-0.1.0.tar.gz

The scripts for the 20/20 Permutation Test can then be found in `Your_Python_Root_Dir/bin`.

Local installation
~~~~~~~~~~~~~~~~~~

Local installation is a good option if you do not have privilege to install a python package and already have the required packages.

**Required packages:**

* numpy
* scipy
* matplotlib (optional, for simulations)
* pandas==0.12.0
* pysam

If you don't have the above required packages, you will need to install them. For the following commands to work you will need `pip <http://pip.readthedocs.org/en/latest/installing.html>`_. If you are using a system wide python, you will need to use `sudo` before the pip command.

.. code-block:: bash

    $ cd permutation2020
    $ pip install -r requirements.txt

Next you will need to build the 20/20 permutation test source files. This is can be accomplished in one command.

.. code-block:: bash

    $ make build

Once finished building you can then use the scripts in the `permutation2020/bin` directory.
