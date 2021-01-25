# Short script to compute the 948 x 948 covariance matrix from the
# RegMap genotype data. In genetics, this is sometimes called the
# "kinship matrix".
library(data.table)

# So that the results are reproducible, set the seed.
set.seed(1)

# Load the genotype data as an n x m matrix, where n is the number of
# genetic markers and m is the number of samples.
cat("Reading genotype data from geno.csv.gz.\n")
geno <- fread("geno.csv.gz",sep = ",",header = TRUE)
geno <- as.matrix(geno)
storage.mode(geno) <- "double"

# Select a random subset of 200,000 genetic markers.
n <- 2e5
i <- sample(ncol(geno),n)
geno <- geno[,i]

# Compute the covariance matrix.
n <- nrow(geno)
m <- ncol(geno)
cat(sprintf("Computing %d x %d covariance matrix using %d markers.\n",n,n,m))
t0   <- proc.time()
geno <- scale(geno,center = TRUE,scale = FALSE)
R    <- tcrossprod(geno)
t1   <- proc.time()
print(t1 - t0)

# Save the covariance matrix.
cat("Saving the covariance matrix to cov.Rdata.\n")
save(list = "R",file = "cov.RData")
