#!/bin/bash
#SBATCH --account=def-stinch     # replace this with your own account
#SBATCH --cpus-per-task=1        # num cpus
#SBATCH --mem-per-cpu=2G         # memory; default unit is megabytes
#SBATCH --time=0-00:15           # time (DD-HH:MM)
module load StdEnv/2020
module load gcc/9.3.0 gurobi/11.0.0 r/4.3.1 

echo "Threads ${SLURM_CPUS_ON_NODE:-1}" > gurobi.env
Rscript wrap_cvxrtest.R