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

The AlphaImpute.zip file contains executibles for Windows, Linux, and Mac, a manual, and an example dataset.

Conditions of use
=================

AlphaImpute is part of a suite of software that our group has developed. It is fully and freely available for all use under the MIT License.

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



Building AlphaImpute
=================

To build alphaimpute please update all of the submodules with `git submodule update --init`. 
Then, use cmake to generate the build system  script with `cmake .`

UNIX-like systems
--------------------

By default, on unix-like systems a make file is generated - this can be built with `make -j8` to build across 8 threads.

Windows systems
--------------------

On windows systems, a visual studio solution will be generated. This can be built in visual studio - be careful!

On windows, 32-bit builds are default - to build for a x64 bit environment, the following has to be done:

`cmake -G "Visual Studio 15 Win64" .`

Also, cmake cannot set the flags correctly for visual studio in fortran - this will have to be done manually.


Build Configurations
--------------------

You can build with debug flags, by using the command `cmake -DCMAKE_BUILD_TYPE=DEBUG .` This allows the program to be debugged with GDB or lldb.

There is also a testing config - which has the fastest performance. This works great when compiling on the same machine, but can cause issues when distributing the binaries. This can be built with `cmake -DCMAKE_BUILD_TYPE=TESTING .` 
