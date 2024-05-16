#!/bin/bash
#SBATCH --account=def-stinch     # replace this with your own account
#SBATCH --nodes=1                # add this line to make sure that slurm uses multiple nodes
#SBATCH --cpus-per-task=36       # number of processes
#SBATCH --mem-per-cpu=500        # memory; default unit is megabytes
#SBATCH --time=0-04:00           # time (DD-HH:MM)
#SBATCH --mail-user=turner.silverthorne@utoronto.ca
#SBATCH --mail-type=ALL
module load StdEnv/2020
module load gcc/9.3.0 r/4.3.1 

cd ~/scratch/oscDetectPaper
Rscript figs_for_paper/Fig1_aug_vs_cvxr/run_calc.R 