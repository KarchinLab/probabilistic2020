Installation
============

Downloading Probabilistic 20/20
-------------------------------

The source files can be cloned from github using `git <http://git-scm.com/>`_:

.. code-block:: bash

    $ git clone https://github.com/ctokheim/probabilistic2020.git

The source files can also be manually downloaded from github at https://github.com/ctokheim/probabilistic2020.

Probabilstic 20/20 can be installed either locally or into your python distribution as a package. 

Python Package Installation
---------------------------

Using the python package installation, all the required python packages for the 20/20 permutation test will automatically be installed for you.

If you are using a system wide python for installation, you will use the following command.

.. code-block:: bash

    $ sudo pip install probabilistic2020-0.1.0.tar.gz

If your python is locally installed then you do not need to use `sudo`.

.. code-block:: bash

    $ pip install probabilistic2020-0.1.0.tar.gz

The scripts for Probabilstic 20/20 can then be found in `Your_Python_Root_Dir/bin`.

Local installation
------------------

Local installation is a good option if you do not have privilege to install a python package and already have the required packages.

**Required packages:**

* numpy
* scipy
* matplotlib (optional, for simulations)
* pandas==0.12.0
* pysam

If you don't have the above required packages, you will need to install them. For the following commands to work you will need `pip <http://pip.readthedocs.org/en/latest/installing.html>`_. If you are using a system wide python, you will need to use `sudo` before the pip command.

.. code-block:: bash

    $ cd probabilstic2020
    $ pip install -r requirements.txt

Next you will need to build the Probabilistic 20/20 source files. This is can be accomplished in one command.

.. code-block:: bash

    $ make build

Once finished building you can then use the scripts in the `probabilstic2020/bin` directory.
