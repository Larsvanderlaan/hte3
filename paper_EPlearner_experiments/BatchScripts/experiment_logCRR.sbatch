#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --output=Rout/par-%J.out
#SBATCH --error=Rout/par-%J.err
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e '
.libPaths( c( "~/Rlibs2") );
n = Sys.getenv("n");
pos = Sys.getenv("pos");
hard = Sys.getenv("hard");
base_learner = Sys.getenv("base_learner");
source("~/EPlearner/FinalSimulationCode/experiment_logCRR.R");
 '
