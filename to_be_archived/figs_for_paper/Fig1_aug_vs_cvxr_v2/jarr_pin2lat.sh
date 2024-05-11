#!/bin/bash
#SBATCH --account=def-stinch   # replace this with your own account
#SBATCH --array=16,24,32       # representing measurement times
#SBATCH --cpus-per-task=10     # number of processes
#SBATCH --mem-per-cpu=150      # memory; default unit is megabytes
#SBATCH --time=0-01:00         # time (DD-HH:MM)
#SBATCH --mail-user=turner.silverthorne@utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load gcc/9.3.0 r/4.3.1  

cd ~/scratch/oscDetectPaper
Rscript figs_for_paper/Fig1_aug_vs_cvxr_v2/sol_pin2lat.R 