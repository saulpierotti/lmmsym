################################################################################
# Example parameter file to be loaded by lmmsym.R
# example: ./lmmsym.R params_file.R
################################################################################

seed <- 1
set.seed(seed)

n_snps_rand <- 1e3

# separate populations
maf_rand <- list(
  runif(n = n_snps_rand, min = 0.2, max = 0.3),
  runif(n = n_snps_rand, min = 0.7, max = 0.8)
)

maf_fixed <- list(
  c(0.2),
  c(0.8)
)

n_samples_l <- list(1500, 500)

# for simplicity just one fixed effect
b <- 0.3 # slope
a <- 2   # intercept

# one beta for each level of the covariate, except the first since dummy encoded
cov1_b <- c(0.3, -0.1)
cov2_b <- c(1, -2)

qcov1_b <- -0.3
qcov1_mean <- -3
qcov1_sd <- 0.2

h2 <- 0.8
var_scale <- 2

ploidy <- 2