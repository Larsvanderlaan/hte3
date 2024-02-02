#!/bin/usr/env bash

nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2500 5000
do
  for pos in "TRUE" "FALSE"
  do
    for hard in "TRUE" "FALSE"
    do
      for base_learner in "xg" "rf" "gam" "earth"
      do
        for sim_type in "CATEhigh", "CATElow"
        do
          sbatch  --export=n=$n,pos=$pos,hard=$hard,base_learner=$base_learner,sim_type=$sim_type ~/LRRsims/BatchScripts/simsCATE.sbatch
        done
      done
    done
  done
done
