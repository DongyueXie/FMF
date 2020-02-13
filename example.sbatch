#!/bin/bash
#SBATCH --partition=broadwl
#SBATCH --mem=8G
#SBATCH --time=2:00:00
# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.5.1
# Run the R script.
echo "run example"
# export OPENBLAS_NUM_THREADS=1
#R CMD BATCH --no-save --no-restore --slave fit_droplet_altsqp.R \
#  ../output/fit_droplet_altsqp.Rout

Rscript GTEx_Splicing.R

