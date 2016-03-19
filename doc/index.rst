.. Probabilistic 20/20 documentation master file, created by
   sphinx-quickstart on Mon Jul 28 13:53:42 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Probabilistic 20/20
===================

:Author: Collin Tokheim
:Contact: ctokheim@jhu.edu
:License: To be Decided
:Source code: `GitHub <https://github.com/ctokheim/probabilistic2020>`_
:Q&A: `Biostars <https://www.biostars.org/t/prob2020/>`_ 

The Probabibilistic 20/20 test identifies genes with signficant oncogene-like and tumor suppressor gene-like mutational patterns for small coding region variants. 
Putative signficant oncogenes are found through evaluating 
missense mutation clustering and *in silico* pathogenicity scores. Often highly clustered missense
mutations are indicative of activating mutations.
While statistically signficant tumor suppressor genes (TSGs) are found by abnormally high number of inactivating mutations.

Probabilistic 20/20 evaluates statistical significance by employing 
monte carlo simulations, which incorporates observed mutation context. Monte carlo
simulations are performed within the same gene and thus avoid building a background
distribution based on other genes.  

The Probabilistic 20/20 test has nice properties since it accounts
for several factors that could effect the significance of driver genes.

* gene length
* mutation context
* gene sequence (e.g. codon bias)

Contents:

.. toctree::
   :maxdepth: 2

   installation
   introduction
   tutorial
   user_doc
   dev_doc
      api
   faq


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

