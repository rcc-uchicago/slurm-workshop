# Short script to compute the 948 x 948 kinship matrix from the RegMap
# genotype data.

# SET UP ENVIRONMENT
# ------------------
library("data.table")

set.seed(1)

# IMPORT GENOTYPE DATA
# --------------------
# Here the genotype data are loaded as a matrix of floating-point numbers.
cat("Reading genotype data.\n")
geno <- fread("geno.csv.gz",sep = ",",header = TRUE)
geno <- as.matrix(geno)
geno <- t(geno)
storage.mode(geno) <- "double"

n <- 40000
m <- 500
i <- sample(nrow(geno),n)
j <- sample(ncol(geno),m)
geno <- geno[i,j]

# COMPUTE KINSHIP MATRIX
# ----------------------
# The kinship matrix is computed by taking the cross-product of the
# genotype matrix.
n <- ncol(geno)
cat(sprintf("Computing %d x %d covariance matrix.\n",n,n))
timing <- system.time(R <- cov(geno))
print(timing)
