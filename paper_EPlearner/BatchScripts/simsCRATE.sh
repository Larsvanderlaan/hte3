#!/bin/usr/env bash


export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2000 3000 4000 5000
do
  for pos in "TRUE" "FALSE"
  do
    for hard in "TRUE" "FALSE"
    do
      for base_learner in "xg" "rf" "gam" "glm"
      do
        for sim_type in "CRATE"
        do
          sbatch  --export=n=$n,pos=$pos,hard=$hard,base_learner=$base_learner,sim_type=$sim_type ~/EPlearner/BatchScripts/simsCRATE.sbatch
        done
      done
    done
  done
done
