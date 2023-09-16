#!/usr/bin/env Rscript

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
  make_option(
    "--seed",
    type = "integer",
    help = "Random seed for reproducibility",
    default = 1
  ),
  make_option(
    "--n_snps_rand",
    type = "integer",
    help = "Number of SNPs modeled as a random effect",
    default = 1e3
  ),
  make_option(
    "--n_snps_fixed",
    type = "integer",
    help = "Number of SNPs modeled as a fixed effect",
    default = 2
  ),
  make_option(
    "--n_populations",
    type = "integer",
    help = "Number of sub-populations with distinct allele frequencies",
    default = 2
  ),
  make_option(
    "--rand_snp_freq",
    type = "character",
    help = reflow_opt(
      "Frequency for SNPs modeled as random effects.",
      "If comma-separated list, frequency in each sub-population."
    ),
    default = "0.2,0.8"
  ),
  make_option(
    "--fixed_snp_freq",
    type = "character",
    help = reflow_opt(
      "Frequency for SNPs modeled as fixed effects.",
      "If comma-separated list, frequency in each sub-population."
    ),
    default = "0.8,0.2"
  ),
  make_option(
    "--n_samples",
    type = "character",
    help = reflow_opt(
      "Number of samples. If comma-separated list,",
      "number of samples in each sub-population."
    ),
    default = "1500,500"
  ),
  make_option(
    "--beta_fixed",
    type = "character",
    help = reflow_opt(
      "Comma-separated list of beta values for the fixed effect SNPs.",
      "Must be of length equal to 'n_snps_fixed'"
    ),
    default = "0.3,-0.3"
  ),
  make_option(
    "--intercept",
    type = "numeric",
    help = "Intercept of the phenotype",
    default = 2
  ),
  make_option(
    "--ploidy",
    type = "integer",
    help = reflow_opt(
      "Ploidy of the organism.",
      "Affects the range of dosages, which will be [0,ploidy]"
    ),
    default = 2
  ),
  make_option(
    "--heritability",
    type = "numeric",
    help = reflow_opt(
      "Narrow-sense heritability of the simulated phenotype.",
      "Must be between 0 and 1."
    ),
    default = 0.8
  ),
  make_option(
    "--variance_scale",
    type = "numeric",
    help = reflow_opt(
      "Scale of the genetic and environmental variances.",
      "The genetic variance is heritability * variance_scale, while the",
      "environmental variance is (1 - heritability) * variance_scale"
    ),
    default = 2
  ),
  make_option(
    "--n_qcov",
    type = "numeric",
    help = "Number of quantitative covariates",
    default = 2
  ),
  make_option(
    "--beta_qcov",
    type = "character",
    help = reflow_opt(
      "Beta for the quantitative covariate.",
      "If comma-separated list, beta for each of the quantitative",
      "covariates. Must be of length 'n_qcov'"
    ),
    default = "2,-3"
  ),
  make_option(
    "--mean_qcov",
    type = "character",
    help = reflow_opt(
      "Mean for the quantitative covariate.",
      "If comma-separated list, mean for each of the quantitative",
      "covariates. Must be of length 'n_qcov'"
    ),
    default = "-2,3"
  ),
  make_option(
    "--sd_qcov",
    type = "character",
    help = reflow_opt(
      "Standard deviation for the quantitative covariate.",
      "If comma-separated list, standard deviation for each of the",
      "quantitative covariates. Must be of length 'n_qcov'"
    ),
    default = "0.2,0.3"
  ),
  make_option(
    "--n_cov",
    type = "numeric",
    help = "Number of covariates",
    default = 1
  ),
  make_option(
    "--n_levels_cov",
    type = "character",
    help = reflow_opt(
      "Number of levels for the covariate.",
      "If comma-separated list, number of levels for each of the",
      "covariates. Must be of length 'n_cov'. Each covariate must have",
      "> 1 level"
    ),
    default = 2
  ),
  make_option(
    "--beta_cov",
    type = "character",
    help = reflow_opt(
      "Beta for each level of the covariate. Comma-separated list with one",
      "value for each level of the covariate - 1 (dummy encoding is used),",
      "with betas for the levels of multiple covariates given sequentially.",
      "Must be of length equal to the sum of (n_levels - 1) for each covariate."
    ),
    default = "3,5",
    meta = "beta_cov1_lev1,beta_cov1_lev2,beta_cov2_lev1,beta_cov2_lev2"
  ),
  make_option(
    "--out",
    type = "character",
    help = "Output name prefix",
    default = format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

banner <- paste(
  "\nlmmsym: a simple simulator for GWAS data with population structure",
  "\nAuthor: Saul Pierotti\n\n",
  sep = "\n"
)

description <- reflow_desc(
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

opt <- parse_args(
  OptionParser(
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
b <- as.numeric(split_vector(opt[["beta_fixed"]]))
a <- opt[["intercept"]]
ploidy <- opt[["ploidy"]]
h2 <- opt[["heritability"]]
var_scale <- opt[["variance_scale"]]
n_qcov <- opt[["n_qcov"]]
beta_qcov <- as.numeric(split_vector(opt[["beta_qcov"]]))
mean_qcov <- as.numeric(split_vector(opt[["mean_qcov"]]))
sd_qcov <- as.numeric(split_vector(opt[["sd_qcov"]]))
n_cov <- opt[["n_cov"]]
n_levels_cov <- as.integer(split_vector(opt[["n_levels_cov"]]))
beta_cov <- as.integer(split_vector(opt[["n_levels_cov"]]))
out <- opt[["out"]]

################################################################################
# Checks
################################################################################

stopifnot(n_populations == length(rand_snp_freq))
stopifnot(n_populations == length(fixed_snp_freq))
stopifnot(n_populations == length(n_samples))
stopifnot(n_snps_fixed == length(b))
stopifnot(n_qcov == length(beta_qcov))
stopifnot(n_qcov == length(mean_qcov))
stopifnot(n_qcov == length(sd_qcov))
stopifnot(n_cov == length(n_levels_cov))
stopifnot(sum(n_levels_cov) - n_cov == length(beta_cov))
stopifnot(sum(n_levels_cov) - n_cov == length(beta_cov))

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
Z.std <- scale(Z)

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

# qcov1 <- matrix(rnorm(n = n_samples, mean = qcov1_mean, sd = qcov1_sd), ncol = 1)
# cov1_f <- as.factor(
#  sample(0:length(cov1_b), size = n_samples, replace = TRUE)
# )
# cov1 <- model.matrix(~cov1_f)[, 2:length(levels(cov1_f))]
# cov2_f <- as.factor(
#  sample(0:length(cov2_b), size = n_samples, replace = TRUE)
# )
# cov2 <- model.matrix(~cov2_f)[, 2:length(levels(cov2_f))]

message("Building phenotype")
y <- a + X %*% b + g + e
stopifnot(dim(y) == c(n_samples_tot, 1))

message("Estimating variance components...")
intercept_col <- matrix(rep(1, n_samples_tot), ncol = 1)
C <- cbind(intercept_col)
gaston_res <- gaston::lmm.aireml(Y = y, X = C, K = K)
s2e_reml <- gaston_res[["sigma2"]]
s2g_reml <- gaston_res[["tau"]]
h2_reml <- s2g_reml / (s2e_reml + s2g_reml)
message("\n*** It is important that the following estimates are accurate ***")
message("h2 realised: ", h2_reml, " true: ", h2)
message("Genetic variance realised: ", s2g_reml, " true: ", var_g)
message("Residual variance realised: ", s2e_reml, " true: ", var_e)
message("***\n")
#
message("Decorrelating")
V <- K * s2g_reml + s2e_reml * diag(1, dim(K))
L <- t(chol(V))
X_mm <- forwardsolve(L, X)
C_mm <- forwardsolve(L, C)
y_mm <- forwardsolve(L, y)


# message("Fitting mixed model")
# fit <- lm(y_mm ~ 0 + C_mm + X_mm)
# a_est <- coef(fit)[[1]]
# b_est <- coef(fit)[[1]]
# message("\n*** It is important that the following estimates are accurate ***")
# message("Intercept estimated: ", a_est, " true: ", a)
# message("Slope estimated: ", b_est, " true: ", b)
# message("\n*** ***")
##
# message("Fitting linear model")
# fit_linear <- lm(y ~ 0 + X + C)
# a_est_linear <- coef(fit_linear)[[2]]
# b_est_linear <- coef(fit_linear)[[1]]
# message(
#  "\n*** The following estimate will probably be off if there ",
#  "is population structure ***"
# )
# message("Intercept estimated: ", a_est_linear, " true: ", a)
# message("Slope estimated: ", b_est_linear, " true: ", b)
#
#################################################################################
## Save output
#################################################################################
#
# res <- list(
#  seed = seed,
#  n_snps_rand = n_snps_rand,
#  n_snps_fixed = n_snps_fixed,
#  maf_rand = maf_rand,
#  maf_fixed = maf_fixed,
#  n_samples_l = n_samples_l,
#  n_samples = n_samples,
#  b = b,
#  a = a,
#  cov1_b - cov1_b,
#  cov2_b = cov2_b,
#  qcov1_b = qcov1_b,
#  qcov1_mean = qcov1_mean,
#  qcov1_sd = qcov1_sd,
#  ploidy = ploidy,
#  h2 = h2,
#  var_scale = var_scale,
#  var_g = var_g,
#  var_e = var_e,
#  K = K,
#  Z = Z,
#  X = X,
#  C = C,
#  Z.std = Z.std,
#  u = u,
#  g = g,
#  y = y,
#  gaston_res = gaston_res,
#  V = V,
#  L = L,
#  X_mm = X_mm,
#  y_mm = y_mm,
#  C_mm = C_mm,
#  fit = fit,
#  a_est = a_est,
#  b_est = b_est,
#  fit_linear = fit_linear,
#  a_est_linear = a_est_linear,
#  b_est_linear = b_est_linear
# )
#
# message("Saving R objects")
# saveRDS(res, sprintf("%s.lmmsym.rds", now))
#
# if (require(ComplexHeatmap) & require(circlize)) {
#  message("Saving relatedness matrix heatmap")
#  png(sprintf("%s.lmmsym.relatedness_matrix.png", now))
#  K.clust <- hclust(as.dist(-K))
#  ComplexHeatmap::Heatmap(
#    unname(K),
#    cluster_rows = K.clust,
#    cluster_columns = K.clust,
#    col = circlize::colorRamp2(c(min(K), mean(K), max(K)), c("blue", "white", "red")),
#    show_row_dend = FALSE,
#    show_column_dend = FALSE,
#    heatmap_legend_param = list(title = "Relatedness"),
#    use_raster = TRUE
#  )
#  dev.off()
# } else {
#  message("ComplexHeatmap and circlize packages not detetcted. Install them if ypu want to plot an heatmap of the relatedness matrix.")
# }
#
# sample_names <- sprintf("sample_%s", 1:n_samples)
#
# message("Saving VCF of the genotypes")
# if (require(vcfR)) {
#  make_fix_line <- function(pos, chr) {
#    c(chr, pos, NA, "A", "T", NA, "PASS", "")
#  }
#  fix_rand <- t(sapply(1:n_snps_rand, make_fix_line, chr = "rand"))
#  fix_fixed <- t(sapply(1:n_snps_fixed, make_fix_line, chr = "fixed"))
#  fix <- rbind(fix_rand, fix_fixed)
#  colnames(fix) <- c(
#    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
#  )
#
#  make_gt <- function(X) {
#    gt <- X
#    gt[X == 0] <- "0/0"
#    gt[X == 1] <- "0/1"
#    gt[X == 2] <- "1/1"
#    return(gt)
#  }
#
#  gt <- t(cbind(make_gt(Z), make_gt(X)))
#  gt <- cbind("GT", gt)
#  colnames(gt) <- c("FORMAT", sample_names)
#
#  vcf_obj <- new(
#    "vcfR",
#    meta = c(
#      "##fileformat=VCFv4.0",
#      "##FILTER=<ID=PASS,Description=\"All filters passed\">",
#      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Hard coded genotype\">",
#      "##contig=<ID=rand>",
#      "##contig=<ID=fixed>",
#      sprintf("##synthetic genetic data created by lmmsym on %s", now)
#    ),
#    fix = fix,
#    gt = gt
#  )
#
#  write.vcf(vcf_obj, sprintf("%s.lmmsym.vcf.gz", now))
# }
#
# message("Saving phenotype file")
# pheno <- data.frame(
#  IID = sample_names,
#  simulated_pheno = y
# )
# write.table(pheno, sprintf("%s.symlmm.pheno", now), row.names = FALSE, quote = FALSE)
#
# message("Saving covariate file")
# covar <- data.frame(
#  IID = sample_names,
#  cov1 = cov1_f,
#  cov2 = cov2_f
# )
# write.table(covar, sprintf("%s.symlmm.covar", now), row.names = FALSE, quote = FALSE)
#
# message("Saving quantitative covariate file")
# qcovar <- data.frame(
#  IID = sample_names,
#  qcov1 = qcov1
# )
# write.table(qcovar, sprintf("%s.symlmm.qcovar", now), row.names = FALSE, quote = FALSE)
