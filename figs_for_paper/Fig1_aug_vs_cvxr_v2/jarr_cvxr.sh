#!/bin/bash
#SBATCH --account=def-stinch     # replace this with your own account
#SBATCH --array=16,24
#SBATCH --mem-per-cpu=8G        # memory; default unit is megabytes
#SBATCH --time=0-00:08 # time (DD-HH:MM)

module load StdEnv/2020
module load gurobi
module load gcc/9.3.0 r/4.3.1  

cd ~/scratch/oscDetectPaper
Rscript figs_for_paper/Fig1_aug_vs_cvxr_v2/sol_cvxr.R 