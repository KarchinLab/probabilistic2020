.. _tutorial-ref:

Tutorial
========

Probabilistic 20/20 consists of two broad statistical tests (oncogene-like and tsg-like) 
and somatic mutation simulation framework. Internally, the simulation framework is 
used to establish statistical significance in the hypothesis test through the 
**probabilistic2020** command. However, the simulation framework through the **mut_annotate** command can 
also be used to create a simulated MAF file where aferwords all mutations are distributed
like passengers based on uniform null distribution. Moreover, a set of mutational
features for each gene representative of driver genes (used in 20/20+) can also be
created.

Input formats
-------------

Mutations
+++++++++

Mutations are provided in a Mutation Annotation Format (MAF) file. 
Columns can be in any order, and only a few columns in the MAF file
are needed. The following is a list of the required columns.

* Hugo_Symbol (or named "Gene")
* Chromosome
* Start_Position
* End_Position
* Reference_Allele
* Tumor_Seq_Allele2 (or named "Tumor_Allele")
* Tumor_Sample_Barcode (or named "Tumor_Sample")
* Variant_Classification

The remaining columns in the MAF specification can be 
left empty or not included. Since a MAF file has many additional 
annotation columns, removing additional columns will reduce
the memory usage of probabilistic2020.

Only coding variants found in the Variant_Classification column will be used, which includes the following: 'Missense_Mutation', 'Silent', 'Nonsense_Mutation', 'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins', 'In_Frame_Del', 'Frame_Shift_Indel', or 'In_Frame_Indel'. Note, although 'In_Frame_Indel' and 'Frame_Shift_Indel' are not official MAF specification values, for the purpose of this program represent either and insertion or deletion. Other values for the Variant_Classification column are assumed to be non-coding, and dropped from the analysis.

Gene BED file
+++++++++++++

A single reference transcript for each gene is stored in BED12 format. Instead of
using the transcript name for the name field in the BED file,
the gene symbol which matches the MAF file should be used.
In the example data, the longest CDS transcript from SNVBox was used.

.. _make-fasta:

Gene FASTA
++++++++++

Gene sequences are extracted from a genome FASTA file, and is a step that only needs to be done once.  
To do this, you need a BED file with names corresponding to genes, and a genome FASTA (e.g. hg19).
You can download hg19 from `here <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit>`_.
Creating the gene sequence FASTA is then done by the `extract_gene_seq` script:

.. code-block:: bash

    $ extract_gene_seq -i hg19.fa -b snvboxGenes.bed -o snvboxGenes.fa

In this case the BED file is created using SNVBox, a genome FASTA file for hg19 (hg19.fa), and the
resulting coding sequences for the gene are stored in snvboxGenes.fa.

Pre-computed scores (optional)
++++++++++++++++++++++++++++++

Two pre-computed scores are used to evaluate missense pathogenicity 
scores and evolutionary conservation. Both are provided in the example
data, matching the reference transcript annotation from SNVBox.
Including the score information is useful, but optional. The 
pre-computed missense pathogenicity scores are from the 
`VEST algorithm <http://www.ncbi.nlm.nih.gov/pubmed/23819870>`_.
The evolutionary conservation scores are calculated as the entropy of 
a specific column in the protein-translated version of UCSC's 46-way vertebrate alignment.

Running the statistical test
----------------------------

The statistical tests account for gene sequence and mutational base context.
Each gene is represented by a single reference transcript (above is longest CDS SNVBox transcript).
By default the relevant sequence context for mutations are utilized from
`CHASM paper <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2763410/>`_ (denoted by **-c 1.5** parameter). This includes some common dinucletoide contexts
like CpG, and otherwise just a single base. Ultimately a multiple testing corrected q-value
is reported using the Benjamini-Hochberg (BH) method.

**Technical detail:** Running on the obtained pan-cancer data may take several hours to run on a single
core. Specifying the **-p** parameter to use multiple processors will speed up run time if available.
Lowering the number of iterations (default: 100,000) will decrease run time, but also decrease the resolution
of p-values.

Running oncogene sub-command
++++++++++++++++++++++++++++

The oncogene sub-command examines missense position clustering (by codon) and elevated
*in silico* pathogenicity scores (VEST). The score directory contains pre-computed values for VEST scores.
The p-values will be combined using fisher's method
to report a single p-value with a BH FDR. In the below example, the command is parallelized
onto 10 processors with the **-p** parameter. Lower this if the compute is not available.

.. code-block:: bash

   $ probabilistic2020 oncogene \
        -i genes.fa \
        -b genes.bed \
        -s score_dir \
        -m mutations.txt \
        -c 1.5
        -p 10 \
        -o oncogene_output.txt

Where genes.fa is your gene FASTA file for your reference transcripts in genes.bed, mutations.txt is your MAF file containing mutations, score_dir is the directory containing the pre-computed VEST scores, and oncogene_output.txt is the file name to save the results.

Output format
#############

The oncogene statistical test will output a tab-delimited file having columns for the 
p-values and Benjamini-Hochberg q-values:

* "entropy"
* "vest" (only included if score_dir provided)
* "combined" (only included if score_dir provided)

The entropy columns evaluate missense clustering at the same codon by using a normalized missense position entropy statistic. Low values for entropy correspond to increased clustering
of missense mutations. The vest columns examine whether missense mutations tend to have
higher *in silico* pathogenicity scores for missense mutations than expected. The "combined"
columns, combine the p-values from VEST scores and missense clustering using Fisher's method.

Running tsg sub-command
+++++++++++++++++++++++

The **tsg** sub-command evaluates for elevated proportion of inactivating point mutations to find TSG-like genes.

.. code-block:: bash

   $ probabilistic2020 tsg \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -p 10 \
        -c 1.5 \
        -o tsg_output.txt

Where genes.fa is your gene FASTA file for your reference transcripts in genes.bed, mutations.txt is your MAF file containing mutations, and tsg_output.txt is the file name to save the results.

Output format
#############

The tsg statistical test examines inactivating single nucleotide variants (nonsense, 
splice site, lost start, and lost stop). Both the p-value ("inactivating p-value")
and the Benjamini-hochberg q-value ("inactivating BH q-value") are reported for 
a higher than expected fraction of inactivating mutations. Mutations which could
not be placed onto the reference transcript will be indicated in the 
"SNVs Unmapped to Ref Tx" column.

Running hotmaps1d sub-command
+++++++++++++++++++++++++++++

The **hotmaps1d** sub-command evaluates particular amino acid residues for elevated cluster of missense mutations in the protein sequence.

.. code-block:: bash

   $ probabilistic2020 hotmaps1d \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -w 3 \
        -p 10 \
        -c 1.5 \
        -o hotmaps1d_output.txt

Where genes.fa is your gene FASTA file for your reference transcripts in genes.bed, mutations.txt is your MAF file containing mutations, and hotmaps1d_output.txt is the file name to save the results. HotMAPS 1D also takes a window size for examining missense mutation clustering. In the above example, the parameter **-w 3** considers 3 residues on either side of each mutated residue. A large number of mutations in this small window may indicate the mutations form a "hotspot", and likely contain driver mutations at the mutated residue. The window size can be changed depending on the preferred granularity of the analysis.

Output format
#############

The hotmaps1d statistical test examines the position of missense mutations in sequence. 
Both the p-value ("p-value") and the Benjamini-hochberg q-value ("q-value") are reported for 
a higher than expected ammount of missense mutations within a given window around a mutation. The "mutation count" column reports how many missense mutations were observed at the particular codon, and the "windowed sum" column reports how many missense mutations were observed in a sequence window encompassing the particular codon.

Simulating somatic mutations
----------------------------

The probabilistic2020 package also allows saving the results of underlying simulation
of somatic mutations. The simulations need a set of observed mutations to create simulated 
mutations. Briefly, for each gene, SNVs (single nucleotide variants) are moved with uniform probability to any matching position in the gene sequence, holding the total number of SNVs fixed.  A matching position was required to have the same base context (e.g. **-c 1.5** = C\*pG, CpG\*, TpC\*, G\*pA, A, C, G, T) as the observed position.  This method of generating a null distribution controls for the particular gene sequence, gene length and mutation base context.  
To simulate small insertions/deletions (indels), indels are moved to different genes according to a multinomial model where the probability is proprotional to the gene length.
This can be done for both creating a simulated MAF file or simulated
features calculated from the mutations.

Simulations are performed with the **mut_annotate** command. The **--seed** parameter
will pass a seed to the pseudo random number generator. If you are performing several
simulations for MAF files and features, then it is critical that every time the seed for each
simulation match. 

Simulated MAF
+++++++++++++

MAF output is designated with the **--maf** flag, but is a substantially reduced version 
then a typical MAF file because it only contains the relevant columns noted in the
mutations input format section. To indicate mutations for each gene should be simulated
once, the **-n 1** parameter is used. If zero is supplied for this parameter, then
simulations are not performed and rather the observed mutations are just annotated
as a MAF file on the corresponding reference transcripts in genes.bed. The pseudo random
number generator seed can be passed with the **--seed** argument.

.. code-block:: bash

   $ mut_annotate \
        --maf \
        -n 1 \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -p 10 \
        -c 1.5 \
        -o maf_output.txt


Simulated Features
++++++++++++++++++

Simulated features which serve as input to `20/20+ <http://2020plus.readthedocs.io/>`_
can also be generated.

.. code-block:: bash

   $ mut_annotate \
        --summary \
        -n 1 \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -p 10 \
        -c 1.5 \
        -o summary_output.txt
