################################################################################
# Input parameters
################################################################################
set.seed(1)

n_snps_rand <- 1e4
maf_rand <- runif(n = n_snps_rand, min = 0.01, max = 0.5)
b <- 3 # beta coefficient
a <- 2 # intercept
ploidy <- 2
n_founders <- 8
cross_structure <- list(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(1, 6),
  c(1, 7),
  c(1, 8),
  c(2, 3),
  c(4, 5),
  c(6, 7)
)
samples_per_cross <- c(100, 100, 500, 500, 100, 100, 100, 100, 100, 300)
n_crossing_over <- 3
h2 <- 0.8
var_scale <- 2

################################################################################
# Compute values
################################################################################

n_crosses <- length(cross_structure)
n_samples <- sum(samples_per_cross)
var_g <- h2 * var_scale
var_e <- (1 - h2) * var_scale

make_chromosome <- function(n_snps, maf) {
  founder_haps <- sapply(
    maf,
    function(p) {
      rbinom(n = n_founders, size = 1, prob = p)
    }
  )
  
  stopifnot(dim(founder_haps) == c(n_founders, n_snps_rand))
  stopifnot(length(cross_structure) == length(samples_per_cross))
  stopifnot(length(cross_structure) == length(unique(cross_structure)))
  
  X <- numeric(0)
  for (cross_idx in 1:n_crosses) {
    cross_founders <- cross_structure[[cross_idx]]
    hap1_idx <- cross_founders[[1]]
    hap2_idx <- cross_founders[[2]]
    hap1 <- founder_haps[hap1_idx, ]
    hap2 <- founder_haps[hap2_idx, ]

    for (sample_idx in 1:samples_per_cross[cross_idx]) {
      co_events <- sort(as.integer(
        round(runif(n = n_crossing_over, min = 2, max = n_snps_rand - 2))
      ))

      # this is the A/B state of the F2
      states1 <- rbinom(n = n_crossing_over + 1, size = 1, prob = 0.5)
      states2 <- rbinom(n = n_crossing_over + 1, size = 1, prob = 0.5)

      f2_dosage <- numeric(0)
      for (co_idx in 1:(n_crossing_over + 1)) {
        start <- if (co_idx > 1) co_events[co_idx - 1] else 1
        end <- if (co_idx <= n_crossing_over) {
          co_events[co_idx] - 1
        } else {
          n_snps_rand
        }
        f2_dosage <- c(
          f2_dosage,
          hap1[start:end] * states1[co_idx] + hap2[start:end] * states2[co_idx]
        )
      }
      message("Generatyng dosage cross: ", cross_idx, ", sample: ", sample_idx)
      stopifnot(length(f2_dosage) == n_snps_rand)
      stopifnot(all(!is.na(f2_dosage)))
      stopifnot(all(f2_dosage %in% 0:2))
      X <- rbind(X, f2_dosage)
    }
  }
  return(X)
}

Z <- make_chromosome(n_snps_rand, maf_rand)

stopifnot(dim(Z) == c(n_samples, n_snps_rand))

Z.std <- (Z - (ploidy * maf_rand)) / sqrt(ploidy * maf_rand * (1 - maf_rand))
K <- Z.std %*% t(Z.std)

stopifnot(dim(K) == c(n_samples, n_samples))

u <- rnorm(n_snps, mean = 0, sd = sqrt(var_g / n_snps))
g <- Z_std %*% u
e <- rnorm(n = n_samples, mean = 0, sd = sqrt(var_e))
y <- a + X %*% b + g + e
