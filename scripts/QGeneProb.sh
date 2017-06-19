#!/bin/sh
########################################
#                                      #
# GE job script for ECDF Cluster       #
#                                      #
# by Roberto Antolin                   #
# AlphaGenes Group                     #
#                                      #
########################################

# Grid Engine options
#$ -cwd
#$ -j y
#$ -l h_rt=00:30:00
#$ -l h_vmem=8G

# Initialise the module environment
source /etc/profile.d/modules.sh

# Load modules
module load intel

# Standard report

echo "Working directory:"
pwd
date

#echo
#echo "System PATH (default):"
#echo $PATH
#echo "System PATH modified:"
export PATH=.:~/bin:$PATH
#echo $PATH
#echo "Ulimit:"
#ulimit -a
#echo


time ./geneProb
