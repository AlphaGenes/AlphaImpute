#!/bin/bash

RS=$1
if [[ $RS == 1 ]]; then
    RESTART=""
    if [[ $2 == "Iterate" ]]; then
        ITERATE=$2
    else
        ITERATE=""
    fi
else
    RESTART=$2
fi

WAIT=""

function runGP {
  Iterate=$1
  if [[ $Iterate == "Iterate" ]]; then
      WAIT=""
      if [[ $2 -ne "" ]]; then
        WAITING="-hold_jid $2"
      else
        WAITING=""
      fi
      PREFIX=IGP
  else
      WAITING=""
      PREFIX=GP
  fi

  for GProbs in ${Iterate}GeneProb/GeneProb*; do
    cd $GProbs
    cp ../../scripts/QGeneProb.sh .
    cp ../../geneProb .
    job=$PREFIX$(echo $GProbs | cut -d'/' -f2 | cut -d'b' -f2)
    WAIT="${WAIT},$job"
    #echo "qsub -N $job ./QGeneProb.sh"
    qsub $WAITING -N $job ./QGeneProb.sh
    cd -
  done
}


if [ $RS == 1 ]; then
    runGP $ITERATE
fi