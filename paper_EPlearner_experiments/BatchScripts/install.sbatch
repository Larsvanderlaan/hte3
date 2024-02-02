#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1

module load R
Rscript -e ' source("~/EPlearner/BatchScripts/install.R")'
