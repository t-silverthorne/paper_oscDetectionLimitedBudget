#!/bin/bash
#SBATCH --account=def-stinch   # replace this with your own account
#SBATCH --array=30,31,32,35,36,37,46,47          # representing measurement times
#SBATCH --cpus-per-task=32     # number of processes
#SBATCH --mem-per-cpu=2G       # memory; default unit is megabytes
#SBATCH --time=0-02:10 # time (DD-HH:MM)
#SBATCH --mail-user=turner.silverthorne@utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020
module load gurobi
module load gcc/9.3.0 r/4.3.1  

source ~/.profile
cd ~/scratch/oscDetectPaper
Rscript figs_for_paper/gurobi_raw/gurobi_raw.R 