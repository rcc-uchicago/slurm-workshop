# Short script to compute the 948 x 948 covariance matrix from the
# RegMap genotype data. In genetics, this is sometimes called the
# "kinship matrix".
#
# This script is intended to be run from the command-line shell. It
# accepts a single argument specifying the number of genetic markers
# to use in the covariance matrix calculation. For example, to compute
# the covariance matrix using a random subset of 1,000 markers, run
# this command:
#
#   Rscript compute_cov.R 1000
#
library(data.table)

# Process the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
m    <- as.numeric(args[1])

# So that the results are reproducible, set the seed.
set.seed(1)

# Load the genotype data as an n x m matrix, where n is the number of
# genetic markers and m is the number of samples.
cat("Reading genotype data from geno.csv.gz.\n")
geno <- fread("geno.csv.gz",sep = ",",header = TRUE)
geno <- as.matrix(geno)
storage.mode(geno) <- "double"

# Select a random subset of m genetic markers.
i    <- sample(ncol(geno),m)
geno <- geno[,i]

# Compute the covariance matrix.
n <- nrow(geno)
cat(sprintf("Computing %d x %d covariance matrix using %d markers.\n",n,n,m))
t0   <- proc.time()
geno <- scale(geno,center = TRUE,scale = FALSE)
R    <- tcrossprod(geno)
t1   <- proc.time()
print(t1 - t0)

# Save the covariance matrix.
outfile <- sprintf("cov_%d.RData",m)
cat("Saving the covariance matrix to ",outfile,".\n",sep="")
save(list = "R",file = outfile)
