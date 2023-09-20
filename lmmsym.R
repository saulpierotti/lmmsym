#!/usr/bin/env Rscript

# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Saul Pierotti

################################################################################
################################################################################
# lmmsym: a simple simulator of GWAS data with population structure
#
# Author: Saul Pierotti
################################################################################
################################################################################

if (!suppressPackageStartupMessages(require("optparse"))) {
  install.packages("optparse", repos = "http://cran.rstudio.com/")
}

if (!suppressPackageStartupMessages(require("gaston"))) {
  install.packages("gaston", repos = "http://cran.rstudio.com/")
}

if (!suppressPackageStartupMessages(require("ComplexHeatmap"))) {
  if (!require("BiocManager")) {
    install.packages("BiocManager", repos = "http://cran.rstudio.com/")
  }
  BiocManager::install("ComplexHeatmap")
}

if (!suppressPackageStartupMessages(require("circlize"))) {
  install.packages("circlize", repos = "http://cran.rstudio.com/")
}

if (!suppressPackageStartupMessages(require("magick"))) {
  install.packages("magick", repos = "http://cran.rstudio.com/")
}

if (!suppressPackageStartupMessages(require("vcfR"))) {
  install.packages("vcfR", repos = "http://cran.rstudio.com/")
}

################################################################################
# Parse arguments
################################################################################

split_vector <- function(v) {
  strsplit(v, ",")[[1]]
}

reflow_opt <- function(...) {
  paste(
    strwrap(
      paste(...),
      width = getOption("width") - 8 * 2
    ),
    collapse = "\n\t\t"
  )
}

reflow_desc <- function(...) {
  paste(strwrap(paste(...), width = getOption("width")), collapse = "\n")
}

option_list <- list(
  optparse::make_option(
    "--seed",
    type = "integer",
    help = "Random seed for reproducibility",
    default = 1
  ),
  optparse::make_option(
    "--n_snps_rand",
    type = "integer",
    help = "Number of SNPs modeled as a random effect",
    default = 1e3
  ),
  optparse::make_option(
    "--n_snps_fixed",
    type = "integer",
    help = "Number of SNPs modeled as a fixed effect",
    default = 5
  ),
  optparse::make_option(
    "--n_populations",
    type = "integer",
    help = "Number of sub-populations with distinct allele frequencies",
    default = 2
  ),
  optparse::make_option(
    "--rand_snp_freq",
    type = "character",
    help = reflow_opt(
      "Frequency for SNPs modeled as random effects.",
      "If comma-separated list, frequency in each sub-population."
    ),
    default = "0.2,0.8"
  ),
  optparse::make_option(
    "--fixed_snp_freq",
    type = "character",
    help = reflow_opt(
      "Frequency for SNPs modeled as fixed effects.",
      "If comma-separated list, frequency in each sub-population."
    ),
    default = "0.8,0.2"
  ),
  optparse::make_option(
    "--n_samples",
    type = "character",
    help = reflow_opt(
      "Number of samples. If comma-separated list,",
      "number of samples in each sub-population."
    ),
    default = "1500,500"
  ),
  optparse::make_option(
    "--beta_fixed",
    type = "character",
    help = reflow_opt(
      "Comma-separated list of beta values for the fixed effect SNPs.",
      "Must be of length equal to 'n_snps_fixed'"
    ),
    default = "0.3,0,0,0,0"
  ),
  optparse::make_option(
    "--intercept",
    type = "numeric",
    help = "Intercept of the phenotype",
    default = 2
  ),
  optparse::make_option(
    "--ploidy",
    type = "integer",
    help = reflow_opt(
      "Ploidy of the organism.",
      "Affects the range of dosages, which will be [0,ploidy]"
    ),
    default = 2
  ),
  optparse::make_option(
    "--heritability",
    type = "numeric",
    help = reflow_opt(
      "Narrow-sense heritability of the simulated phenotype.",
      "Must be between 0 and 1."
    ),
    default = 0.8
  ),
  optparse::make_option(
    "--variance_scale",
    type = "numeric",
    help = reflow_opt(
      "Scale of the genetic and environmental variances.",
      "The genetic variance is heritability * variance_scale, while the",
      "environmental variance is (1 - heritability) * variance_scale"
    ),
    default = 2
  ),
  optparse::make_option(
    "--n_qcov",
    type = "numeric",
    help = "Number of quantitative covariates",
    default = 2
  ),
  optparse::make_option(
    "--beta_qcov",
    type = "character",
    help = reflow_opt(
      "Beta for the quantitative covariate.",
      "If comma-separated list, beta for each of the quantitative",
      "covariates. Must be of length 'n_qcov'"
    ),
    default = "2,-3"
  ),
  optparse::make_option(
    "--mean_qcov",
    type = "character",
    help = reflow_opt(
      "Mean for the quantitative covariate.",
      "If comma-separated list, mean for each of the quantitative",
      "covariates. Must be of length 'n_qcov'"
    ),
    default = "-2,3"
  ),
  optparse::make_option(
    "--sd_qcov",
    type = "character",
    help = reflow_opt(
      "Standard deviation for the quantitative covariate.",
      "If comma-separated list, standard deviation for each of the",
      "quantitative covariates. Must be of length 'n_qcov'"
    ),
    default = "0.2,0.3"
  ),
  optparse::make_option(
    "--n_cov",
    type = "numeric",
    help = "Number of covariates",
    default = 2
  ),
  optparse::make_option(
    "--n_levels_cov",
    type = "character",
    help = reflow_opt(
      "Number of levels for the covariate.",
      "If comma-separated list, number of levels for each of the",
      "covariates. Must be of length 'n_cov'. Each covariate must have",
      "> 1 level"
    ),
    default = "2,3"
  ),
  optparse::make_option(
    "--beta_cov",
    type = "character",
    help = reflow_opt(
      "Beta for each level of the covariate. Comma-separated list with one",
      "value for each level of the covariate - 1 (dummy encoding is used),",
      "with betas for the levels of multiple covariates given sequentially.",
      "Must be of length equal to the sum of (n_levels - 1) for each covariate."
    ),
    default = "3,5,2",
    meta = "beta_cov1_lev1,beta_cov1_lev2,beta_cov2_lev1,beta_cov2_lev2"
  ),
  optparse::make_option(
    "--out",
    type = "character",
    help = "Output name prefix",
    default = format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

banner <- paste(
  "\nlmmsym v0.1.1: a simple simulator for GWAS data with population structure",
  "Author: Saul Pierotti\n",
  sep = "\n"
)

description <- paste0(
  "\n",
  reflow_desc(
    "A simple simulator for gwas data. It simulates a population",
    "composed of sub-populations, which differ in terms of allele",
    "frequencies. Two chromosomes are created, one called 'rand'",
    "and one called 'fixed'. The idea is that 'rand' contains SNPs",
    "modeled as random effects while 'fixed' SNPs modeled as fixed",
    "effects. If you want to test a scenario where a mixed model works",
    "well but a linear model fails, try to specify very different allele",
    "frequencies for the random and fixed effects in the sub-populations.",
    "Note that some combinations may yield data which is very hard to",
    "model, so always check that the modelling done by the script itself",
    "is accurate before using the synthetic data to test other tools."
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(
    option_list = option_list,
    description = paste0(banner, description)
  )
)

################################################################################
# Initialise variables
################################################################################

set.seed(opt[["seed"]])

n_snps_rand <- opt[["n_snps_rand"]]
n_snps_fixed <- opt[["n_snps_fixed"]]
n_populations <- opt[["n_populations"]]
rand_snp_freq <- as.numeric(split_vector(opt[["rand_snp_freq"]]))
fixed_snp_freq <- as.numeric(split_vector(opt[["fixed_snp_freq"]]))
n_samples <- as.integer(split_vector(opt[["n_samples"]]))
beta_fixed <- as.numeric(split_vector(opt[["beta_fixed"]]))
intercept <- opt[["intercept"]]
ploidy <- opt[["ploidy"]]
h2 <- opt[["heritability"]]
var_scale <- opt[["variance_scale"]]
n_qcov <- opt[["n_qcov"]]
beta_qcov <- as.numeric(split_vector(opt[["beta_qcov"]]))
mean_qcov <- as.numeric(split_vector(opt[["mean_qcov"]]))
sd_qcov <- as.numeric(split_vector(opt[["sd_qcov"]]))
n_cov <- opt[["n_cov"]]
n_levels_cov <- as.integer(split_vector(opt[["n_levels_cov"]]))
beta_cov <- as.integer(split_vector(opt[["beta_cov"]]))
out <- opt[["out"]]

log <- file(sprintf("%s.lmmsym.log", out), "w")
sink(file = log, type = "message")

message(banner)
message("*** Parameters ***")
for (param in names(opt)) {
  message(param, " = ", opt[[param]])
}
message("***\n")

################################################################################
# Checks
################################################################################

message("Validating parameters")
stopifnot(n_populations == length(rand_snp_freq))
stopifnot(n_populations == length(fixed_snp_freq))
stopifnot(n_populations == length(n_samples))
stopifnot(n_snps_fixed == length(beta_fixed))
stopifnot(n_qcov == length(beta_qcov))
stopifnot(n_qcov == length(mean_qcov))
stopifnot(n_qcov == length(sd_qcov))
stopifnot(n_cov == length(n_levels_cov))
stopifnot(sum(n_levels_cov) - n_cov == length(beta_cov))
stopifnot(sum(n_levels_cov) - n_cov == length(beta_cov))
stopifnot(all(n_levels_cov > 1))

#################################################################################
## Compute values
#################################################################################

n_samples_tot <- sum(n_samples)
var_g <- h2 * var_scale
var_e <- (1 - h2) * var_scale

make_chromosome <- function(p, n) {
  ret <- numeric(0)
  for (pop_idx in 1:n_populations) {
    curr_p <- p[[pop_idx]]
    curr_n_samples <- n_samples[[pop_idx]]
    mat <- replicate(
      expr = rbinom(size = ploidy, prob = curr_p, n = curr_n_samples),
      n = n
    )
    ret <- rbind(ret, mat)
  }
  return(ret)
}


message("Generating random effects SNPs")
Z <- make_chromosome(rand_snp_freq, n_snps_rand)
stopifnot(dim(Z) == c(n_samples_tot, n_snps_rand))

message("Generating fixed effects SNPs")
X <- make_chromosome(fixed_snp_freq, n_snps_fixed)
stopifnot(dim(X) == c(n_samples_tot, n_snps_fixed))

message("Scaling random effect SNPs")
Z.p <- colMeans(Z)
Z.std <- t((t(Z) - (2 * Z.p)) / sqrt(2 * Z.p * (1 - Z.p)))

message("Computing relatedness matrix")
K <- (Z.std %*% t(Z.std)) / n_snps_rand
stopifnot(dim(K) == c(n_samples_tot, n_samples_tot))

message("Building random effect")
u <- rnorm(n_snps_rand, mean = 0, sd = sqrt(var_g / n_snps_rand))
g <- Z.std %*% u
stopifnot(dim(g) == c(n_samples_tot, 1))

message("Building error vector")
e <- matrix(rnorm(n = n_samples_tot, mean = 0, sd = sqrt(var_e)), ncol = 1)
stopifnot(dim(e) == c(n_samples_tot, 1))

message("Building covariates")
qcov_m <- sapply(
  1:n_qcov,
  function(i) {
    rnorm(n = n_samples_tot, mean = mean_qcov[[i]], sd = sd_qcov[[i]])
  }
)
stopifnot(dim(qcov_m) == c(n_samples_tot, n_qcov))

cov_l <- lapply(
  1:n_cov,
  function(i) {
    n_lev <- n_levels_cov[[i]]
    m <- t(rmultinom(n = n_samples_tot, size = 1, prob = rep(1 / n_lev, n_lev)))
    # drop first level since modeled from intercept (dummy encoding)
    m[, 2:n_lev]
  }
)
cov_m <- do.call("cbind", cov_l)
stopifnot(dim(cov_m) == c(n_samples_tot, sum(n_levels_cov) - n_cov))

message("Building final design matrix")
intercept_col <- rep(1, n_samples_tot)
X <- cbind(intercept_col, X, qcov_m, cov_m)
snp_names <- sapply(
  1:n_snps_fixed,
  function(i) {
    sprintf("snp%s", i)
  }
)
qcov_names <- sapply(
  1:n_qcov,
  function(i) {
    sprintf("qcov%s", i)
  }
)
cov_names_l <- lapply(
  1:n_cov,
  function(i) {
    sapply(
      2:n_levels_cov[[i]],
      function(l_i, i) {
        sprintf("cov%sl%s", i, l_i)
      },
      i = i
    )
  }
)
cov_names <- unlist(cov_names_l)
colnames(X) <- c("intercept", snp_names, qcov_names, cov_names)
stopifnot(
  dim(X) == c(
    n_samples_tot,
    1 + n_snps_fixed + n_qcov + sum(n_levels_cov) - n_cov
  )
)

# null design matrix for variance estimation
X_null <- cbind(intercept_col, qcov_m, cov_m)
stopifnot(
  dim(X_null) == c(
    n_samples_tot,
    1 + n_qcov + sum(n_levels_cov) - n_cov
  )
)

message("Building fixed effects vector")
b <- c(intercept, beta_fixed, beta_qcov, beta_cov)
stopifnot(length(b) == 1 + n_snps_fixed + n_qcov + sum(n_levels_cov) - n_cov)

message("Building phenotype")
y <- X %*% b + g + e
stopifnot(dim(y) == c(n_samples_tot, 1))

message("Estimating variance components...")
# here the right thing to do would be to set X = X, but I use X = X_null to
# mimic what would be done in a normal GWAS, where the heritability is
# estimated without taking the fixed-effect SNPs into account
gaston_res <- gaston::lmm.aireml(Y = y, X = X_null, K = K, verbose = FALSE)
s2e_reml <- gaston_res[["sigma2"]]
s2g_reml <- gaston_res[["tau"]]
h2_reml <- s2g_reml / (s2e_reml + s2g_reml)
message("\n*** It is important that the following estimates are accurate ***")
message("h2 realised: ", round(h2_reml, 2), " true: ", h2)
message("Genetic variance realised: ", round(s2g_reml, 2), " true: ", var_g)
message("Residual variance realised: ", round(s2e_reml, 2), " true: ", var_e)
message("***\n")
#
message("Decorrelating")
V <- K * s2g_reml + s2e_reml * diag(1, dim(K))
L <- t(chol(V))
X_mm <- forwardsolve(L, X)
y_mm <- forwardsolve(L, y)
colnames(X_mm) <- colnames(X)

message("Fitting mixed model")
fit <- lm(y_mm ~ 0 + X_mm)

intercept_est <- coef(fit)[["X_mmintercept"]]
beta_fixed_est <- coef(fit)[sprintf("X_mm%s", snp_names)]
message("\n*** It is important that the following estimates are accurate ***")
message("Mixed model")
message("Intercept estimated: ", intercept_est, " true: ", intercept)
for (i in seq_along(snp_names)) {
  message(
    "Slope estimated SNP ", i, ": ",
    beta_fixed_est[[i]], " true: ", beta_fixed[[i]]
  )
}
message("***\n")

message("Fitting linear model")
fit <- lm(y ~ 0 + X)

intercept_est <- coef(fit)[["Xintercept"]]
beta_fixed_est <- coef(fit)[sprintf("X%s", snp_names)]
message("\n*** These estimates are expected to be inaccurate ***")
message("Linear model")
message("Intercept estimated: ", intercept_est, " true: ", intercept)
for (i in seq_along(snp_names)) {
  message(
    "Slope estimated SNP ", i, ": ",
    beta_fixed_est[[i]], " true: ", beta_fixed[[i]]
  )
}
message("***\n")

###############################################################################
# Save output
###############################################################################

message("Saving relatedness matrix heatmap")
ret <- png(sprintf("%s.lmmsym.relatedness_matrix.png", out))
K.clust <- hclust(as.dist(-K))
ComplexHeatmap::Heatmap(
  unname(K),
  cluster_rows = K.clust,
  cluster_columns = K.clust,
  col = circlize::colorRamp2(
    c(min(K), mean(K), max(K)), c("blue", "white", "red")
  ),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  heatmap_legend_param = list(title = "Relatedness"),
  use_raster = TRUE
)
ret <- dev.off()
saveRDS(K, sprintf("%s.lmmsym.relmat.rds", out))


message("Saving VCF of the genotypes")
sample_names <- sprintf("sample_%s", 1:n_samples_tot)

make_fix_line <- function(pos, chr) {
  c(chr, pos, NA, "A", "T", NA, "PASS", "")
}
fix_line_rand <- t(sapply(1:n_snps_rand, make_fix_line, chr = "rand"))
fix_line_fixed <- t(sapply(1:n_snps_fixed, make_fix_line, chr = "fixed"))
fix <- rbind(fix_line_rand, fix_line_fixed)
colnames(fix) <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
)

make_gt <- function(dosage_mat) {
  gt <- matrix(
    data = character(1), ncol = ncol(dosage_mat), nrow = nrow(dosage_mat)
  )
  gt[dosage_mat == 0] <- "0/0"
  gt[dosage_mat == 1] <- "0/1"
  gt[dosage_mat == 2] <- "1/1"
  return(gt)
}
gt <- t(cbind(make_gt(Z), make_gt(X)))
gt <- cbind("GT", gt)
colnames(gt) <- c("FORMAT", sample_names)

vcf_obj <- new(
  vcfR::.__C__vcfR,
  meta = c(
    "##fileformat=VCFv4.0",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Hard coded genotype\">",
    "##contig=<ID=rand>",
    "##contig=<ID=fixed>",
    sprintf("##synthetic genetic data created by lmmsym on %s", out)
  ),
  fix = fix,
  gt = gt
)
vcfR::write.vcf(vcf_obj, sprintf("%s.lmmsym.vcf.gz", out))

message("Saving phenotype file")
pheno <- data.frame(
  IID = sample_names,
  simulated_pheno = y
)
write.table(
  pheno, sprintf("%s.lmmsym.pheno", out),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

message("Saving covariate file")
get_factor <- function(v) {
  ifelse(any(v == 1), which(v == 1), 0)
}
start <- 1
cov_f <- numeric(0)
for (i in 1:n_cov) {
  n_dummies <- n_levels_cov[[i]] - 1
  end <- start + n_dummies - 1
  mat <- matrix(cov_m[, start:end], ncol = n_dummies)
  f_vec <- apply(mat, MARGIN = 1, get_factor)
  cov_f <- cbind(cov_f, f_vec)
  start <- end + 1
}
colnames(cov_f) <- sprintf("cov%s", 1:n_cov)
covar <- data.frame(
  IID = sample_names, cov_f
)
write.table(
  covar, sprintf("%s.lmmsym.covar", out),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

message("Saving quantitative covariate file")
colnames(qcov_m) <- sprintf("qcov%s", 1:n_qcov)
qcovar <- data.frame(IID = sample_names, qcov_m)
write.table(
  qcovar, sprintf("%s.lmmsym.qcovar", out),
  row.names = FALSE, quote = FALSE, sep = "\t"
)
