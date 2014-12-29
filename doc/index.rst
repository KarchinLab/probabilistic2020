.. 20/20 Permutation Test documentation master file, created by
   sphinx-quickstart on Mon Jul 28 13:53:42 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to 20/20 Permutation Test's documentation!
==================================================

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

