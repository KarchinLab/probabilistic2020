.. Probabilistic 20/20 documentation master file, created by
   sphinx-quickstart on Mon Jul 28 13:53:42 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Probabilistic 20/20
===================

The Probabibilistic 20/20 test identifies genes with signficant oncogene-like and tumor suppressor gene-like mutational patterns. 
Putative signficant oncogenes are found through evaluating the position of 
missense mutations (clustered missense mutations tend to indicate actiavting mutations).
While statistically signficant tumor suppressor genes (TSGs) are found by abnormally high number of inactivating mutations.

Probabilistic 20/20 evaluates statistical significance by employing 
monte carlo simulations, which incorporates observed mutation context. Monte carlo,
simulations are performed within the same gene and thus avoid building a background
distribution based on other genes.  

The Probabilistic 20/20 test has nice properties since it accounts
for several factors that could effect the significance of driver genes.

* gene length
* mutation context
* codon bias

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

