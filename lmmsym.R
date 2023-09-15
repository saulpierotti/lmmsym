#!/usr/bin/env Rscript

################################################################################
################################################################################
# lmmsym
#
# Simple simulator of GWAS data with population structure
#
# Note: make sure that the LMM fitting in this script is accurate before using
# the output to test other tools! Certain parameter combinations will produce
# very bad estimates (for example, too few samples).
# 
# Author: Saul Pierotti
################################################################################
################################################################################

message("***")
message("This is lmmsym, a simple simulator of GWAS data with population structure")
message("Author: Saul Pierotti")
message("***\n")

################################################################################
# Default parameters
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

################################################################################
# Load custom parameters
################################################################################

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 1){
  message("Loading parameter file: ", args)
  source(args)
} else if (length(args) == 0) {
  message("No arguments provided, proceding with defaults")
} else {
  stop("Only 1 or 0 arguments are allowed. ", length(args), " arguments supplied.")
}

################################################################################
# Compute values
################################################################################

now <- format(Sys.time(), "%Y%m%d_%H%M%S")

n_samples <- sum(as.numeric(n_samples_l))
n_snps_fixed <- length(b)

var_g <- h2 * var_scale
var_e <- (1 - h2) * var_scale

make_chromosome <- function(maf_l) {
  ret <- numeric(0)
  for (pop_idx in 1:length(maf_l)) {
    maf <- maf_l[[pop_idx]]
    curr_n_samples <- n_samples_l[[pop_idx]]
    X <- sapply(
      maf,
      function(curr_maf){
        rbinom(size = ploidy, prob = curr_maf, n = curr_n_samples)
      }
    )
    ret <- rbind(ret, X)
  }
  return(ret)
}

message("Generating random effects SNPs")
Z <- make_chromosome(maf_rand)
stopifnot(dim(Z) == c(n_samples, n_snps_rand))

message("Generating fixed effects SNPs")
X <- make_chromosome(maf_fixed)
stopifnot(dim(X) == c(n_samples, n_snps_fixed))

message("Computing relatedness matrix")
Z.std <- scale(Z)
K <- (Z.std %*% t(Z.std)) / n_snps_rand
stopifnot(dim(K) == c(n_samples, n_samples))

message("Building phenotype")
u <- rnorm(n_snps_rand, mean = 0, sd = sqrt(var_g/n_snps_rand))
g <- Z.std %*% u
e <- rnorm(n = n_samples, mean = 0, sd = sqrt(var_e))
qcov1 <- matrix(rnorm(n = n_samples, mean = qcov1_mean, sd = qcov1_sd), ncol = 1)
cov1_f  <- as.factor(
  sample(0:length(cov1_b), size = n_samples, replace = TRUE)
)
cov1 <- model.matrix(~ cov1_f)[,2:length(levels(cov1_f))]
cov2_f  <- as.factor(
  sample(0:length(cov2_b), size = n_samples, replace = TRUE)
)
cov2 <- model.matrix(~ cov2_f)[,2:length(levels(cov2_f))]
y <- a + X %*% b + qcov1 %*% qcov1_b + cov1 %*% cov1_b + cov2 %*% cov2_b + g + e

intercept <- matrix(rep(1, n_samples), ncol = 1)
C <- cbind(intercept, qcov1, cov1, cov2)

message("Estimating heritability")
gaston_res <- gaston::lmm.aireml(Y = y, X = C, K = K)
s2e_reml <- gaston_res$sigma2
s2g_reml <- gaston_res$tau
h2_reml <- s2g_reml /(s2e_reml + s2g_reml)
message("h2 realised: ", h2_reml, " true: ", h2)
message("Genetic variance realised: ", s2g_reml, " true: ", var_g)
message("Residual variance realised: ", s2e_reml, " true: ", var_e)

message("Decorrelating")
V <- K * s2g_reml + s2e_reml * diag(1, dim(K))
L <- t(chol(V))
X_mm <- forwardsolve(L, X)
C_mm <- forwardsolve(L, C)
y_mm <- forwardsolve(L, y)

message("Fitting LMM")
fit <- lm(y_mm ~ 0 + X_mm + C_mm)
a_est <- coef(fit)[[2]]
b_est <- coef(fit)[[1]]
message("Intercept estimated: ", a_est, " true: ", a)
message("Slope estimated: ", b_est, " true: ", b)

message("Fitting LM")
fit_linear <- lm(y ~ 0 +  X + C)
a_est_linear <- coef(fit_linear)[[2]]
b_est_linear <- coef(fit_linear)[[1]]
message("Intercept estimated: ", a_est_linear, " true: ", a)
message("Slope estimated: ", b_est_linear, " true: ", b)

################################################################################
# Save output
################################################################################

res <- list(
  seed = seed,
  n_snps_rand = n_snps_rand,
  n_snps_fixed = n_snps_fixed,
  maf_rand = maf_rand,
  maf_fixed = maf_fixed,
  n_samples_l = n_samples_l,
  n_samples = n_samples,
  b = b,
  a = a,
  cov1_b - cov1_b,
  cov2_b = cov2_b,
  qcov1_b = qcov1_b,
  qcov1_mean = qcov1_mean,
  qcov1_sd = qcov1_sd,
  ploidy = ploidy,
  h2 = h2,
  var_scale = var_scale,
  var_g = var_g,
  var_e = var_e,
  K = K,
  Z = Z,
  X = X,
  C = C,
  Z.std = Z.std,
  u = u,
  g = g,
  y = y,
  gaston_res = gaston_res,
  V = V,
  L = L,
  X_mm = X_mm,
  y_mm = y_mm,
  C_mm = C_mm,
  fit = fit,
  a_est = a_est,
  b_est = b_est,
  fit_linear = fit_linear,
  a_est_linear = a_est_linear,
  b_est_linear = b_est_linear
)

message("Saving R objects")
saveRDS(res, sprintf("%s.lmmsym.rds", now))

if (require(ComplexHeatmap) & require(circlize)) {
  message("Saving relatedness matrix heatmap")
  png(sprintf("%s.lmmsym.relatedness_matrix.png", now))
  K.clust <- hclust(as.dist(-K))
  ComplexHeatmap::Heatmap(
    unname(K),
    cluster_rows = K.clust,
    cluster_columns = K.clust,
    col = circlize::colorRamp2(c(min(K), mean(K), max(K)), c("blue", "white", "red")),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(title = "Relatedness"),
    use_raster = TRUE
  )
  dev.off()
} else {
  message("ComplexHeatmap and circlize packages not detetcted. Install them if ypu want to plot an heatmap of the relatedness matrix.")
}

sample_names <- sprintf("sample_%s", 1:n_samples)

message("Saving VCF of the genotypes")
if (require(vcfR)) {
  make_fix_line <- function(pos, chr){c(chr, pos, NA, "A", "T", NA, "PASS", "")}
  fix_rand <- t(sapply(1:n_snps_rand, make_fix_line, chr = "rand"))
  fix_fixed <- t(sapply(1:n_snps_fixed, make_fix_line, chr = "fixed"))
  fix <- rbind(fix_rand, fix_fixed)
  colnames(fix) <- c(
    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
  )
  
  make_gt <- function(X){
    gt <- X
    gt[X == 0] <- "0/0"
    gt[X == 1] <- "0/1"
    gt[X == 2] <- "1/1"
    return(gt)
  }
  
  gt <- t(cbind(make_gt(Z), make_gt(X)))
  gt <- cbind("GT", gt)
  colnames(gt) <- c("FORMAT", sample_names)
  
  vcf_obj <- new(
    "vcfR",
    meta = c(
      "##fileformat=VCFv4.0",
      "##FILTER=<ID=PASS,Description=\"All filters passed\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Hard coded genotype\">",
      "##contig=<ID=rand>",
      "##contig=<ID=fixed>",
      sprintf("##synthetic genetic data created by lmmsym on %s", now)
    ),
    fix = fix,
    gt = gt
  )
  
  write.vcf(vcf_obj, sprintf("%s.lmmsym.vcf.gz", now))
}

message("Saving phenotype file")
pheno <- data.frame(
  IID = sample_names,
  simulated_pheno = y
)
write.table(pheno, sprintf("%s.symlmm.pheno", now), row.names = FALSE, quote = FALSE)

message("Saving covariate file")
covar <- data.frame(
  IID = sample_names,
  cov1 = cov1_f,
  cov2 = cov2_f
)
write.table(covar, sprintf("%s.symlmm.covar", now), row.names = FALSE, quote = FALSE)

message("Saving quantitative covariate file")
qcovar <- data.frame(
  IID = sample_names,
  qcov1 = qcov1 
)
write.table(qcovar, sprintf("%s.symlmm.qcovar", now), row.names = FALSE, quote = FALSE)