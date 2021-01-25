# Short script to compute the 948 x 948 covariance matrix from the
# RegMap genotype data. In genetics, this is sometimes called the
# "kinship matrix".
library("data.table")

# So that the results are reproducible, set the seed.
set.seed(1)

# Load the genotype data as an n x m matrix, where n is the number of
# genetic markers and m is the number of samples.
cat("Reading genotype data from geno.csv.gz.\n")
geno <- fread("geno.csv.gz",sep = ",",header = TRUE)
geno <- as.matrix(geno)
geno <- t(geno)
storage.mode(geno) <- "double"

# Select a random subset of 20,000 genetic markers.
n <- 20000
i <- sample(nrow(geno),n)
geno <- geno[i,]

# Compute the covariance matrix.
n <- nrow(geno)
m <- ncol(geno)
cat(sprintf("Computing %d x %d covariance matrix using %d markers.\n",m,m,n))
timing <- system.time(R <- cov(geno))
print(timing)

# Save the covariance matrix.
cat("Saving the covariance matrix to cov.Rdata.\n")
save(list = "R",file = "cov.RData")
