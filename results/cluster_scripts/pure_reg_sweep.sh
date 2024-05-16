#!/bin/bash
#SBATCH --account=def-stinch   # replace this with your own account
#SBATCH --cpus-per-task=32      # number of processes
#SBATCH --mem-per-cpu=1G       # memory; default unit is megabytes
#SBATCH --time=0-00:30         # time (DD-HH:MM)
#SBATCH --mail-user=turner.silverthorne@utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load gurobi
module load gcc/9.3.0 r/4.3.1  

source ~/.profile
cd ~/scratch/oscDetectPaper
Rscript results/cluster_scripts/pure_reg_sweep.R 