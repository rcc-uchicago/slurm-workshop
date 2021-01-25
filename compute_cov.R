# Short script to compute the 948 x 948 kinship matrix from the RegMap
# genotype data.
library("data.table")

set.seed(1)

# Here the genotype data are loaded as a matrix of floating-point numbers.
cat("Reading genotype data.\n")
geno <- fread("geno.csv.gz",sep = ",",header = TRUE)
geno <- as.matrix(geno)
geno <- t(geno)
storage.mode(geno) <- "double"

n <- 20000
i <- sample(nrow(geno),n)
geno <- geno[i,]

# The kinship matrix is computed by taking the cross-product of the
# genotype matrix.
n <- ncol(geno)
cat(sprintf("Computing %d x %d covariance matrix.\n",n,n))
timing <- system.time(R <- cov(geno))
print(timing)
