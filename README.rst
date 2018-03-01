===========
AlphaImpute
===========

AlphaImpute is a software package for imputing and phasing genotype data in populations with pedigree information
available. The program uses segregation analysis and haplotype library imputation (SAHLI) to impute alleles and
genotypes. A complete description of the methods is given in Hickey et al. (2011b). AlphaImpute consists of a single
program however it calls both AlphaPhase1.1 (Hickey et al., 2011a) and GeneProbForAlphaImpute (Kerr and Kinghorn, 1996).
All information on the model of analysis, input files and their layout, is specified in a single parameter file.

Please report bugs or suggestions on how the program / user interface / manual could be improved or made more user
friendly to `John Hickey <John.Hickey@roslin.ed.ac.uk>`_.

Availability
============

AlphaImpute is available from:

https://bitbucket.org/hickeyjohnteam/alphaimpute

Material available comprises the compiled programs for 64 bit Linux and Mac OSX machines, together with a User Manual
and a suite of worked examples.

Conditions of use
=================

AlphaImpute is available to the scientific community free of charge. Users are required, however, to credit its use in
any publications. Commercial users should contact John Hickey.

Suggested Citation
------------------

Hickey et al. (2011). Segregation analysis and haplotype library imputation to impute SNP alleles in pedigreed
populations. Genetics Selection Evolution, 44:9. doi:10.1186/1297-9686-44-9

Disclaimer
==========

While every effort has been made to ensure that AlphaImpute does what it claims to do, there is absolutely no guarantee
that the results provided are correct. Use of AlphaImpute is entirely at your own risk!

Advertisement
=============
Your welcome to check out our Gibbs sampler (AlphaBayes) specifically designed for GWAS and Genomic Selection.
http://sites.google.com/site/hickeyjohn/alphabayes.

Description of methods
======================
The method implemented in AlphaImpute is described in detail in Hickey et al. (2011).

Using AlphaImpute
=================

.. warning:: AlphaImpute only works for multiple chromosomes if using plink files as input.

