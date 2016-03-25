Tutorial
========

Probabilistic 20/20 consists of two broad statistical tests (oncogene-like and tsg-like) and somatic mutation simulation framework. Internally, the simulation framework is 
used to establish statistical significance in the hypothesis test through the 
**probabilistic2020** command. However, the simulation framework through the **mut_annotate** command can 
also be used to create a simulated MAF file where aferwords all mutations are distributed
like passengers based on uniform null distribution. Moreover, a set of mutational
features for each gene representative of driver genes (used in 20/20+) can also be
created.

In this tutorial we go through a somewhat computationally intensive Pan-cancer
analysis. Please see the quick start for a fast example.
