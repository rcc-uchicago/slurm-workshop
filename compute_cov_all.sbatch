#!/bin/bash

#SBATCH --job-name=compute_cov_all
#SBATCH --partition=broadwl
#SBATCH --reservation=slurm_workshop
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# For reproducibility, add version numbers to the "module load" calls.
module load R/3.6.1
export OPENBLAS_NUM_THREADS=4

# Run the R script.
echo "Running compute_cov_all.R"
Rscript compute_cov_all.R > compute_cov_all.Rout
