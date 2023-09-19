# lmmsym

This is a simple R script that simulates genetic dosage and phenotypes. It is possible to specify an arbitrary number of sub-populations, each with a dedicate allele frequency.
Two chromosomes are simulated, one whose SNPs are modelled as fixed effects and one whose SNPs are modelled as random effects.
The intended use of this script is to test linear mixed models implementations for genetics.
If the linear mixed model works properly, it should be able to correctly reconstruct the effect sizes of the fixed effect SNPs and the simulated random effect heritability.
The output is produced in the form of a vcf file for the genotype, and plink-style files for the phenotype and covariates.

## Usage

Navigate to the root of the repository, make the script executable and just run it

```
chmod +x lmmsym.R
./lmmsym.R
```

If you want to customise some parameters of the simulation, explore the available options with

```
./lmmsym.R --help
```

## Parameters

```
Usage: ./lmmsym.R [options]

lmmsym: a simple simulator for GWAS data with population structure
Author: Saul Pierotti

A simple simulator for gwas data. It simulates a population composed of
sub-populations, which differ in terms of allele frequencies. Two chromosomes
are created, one called 'rand' and one called 'fixed'. The idea is that 'rand'
contains SNPs modeled as random effects while 'fixed' SNPs modeled as fixed
effects. If you want to test a scenario where a mixed model works well but a
linear model fails, try to specify very different allele frequencies for the
random and fixed effects in the sub-populations. Note that some combinations
may yield data which is very hard to model, so always check that the modelling
done by the script itself is accurate before using the synthetic data to test
other tools.

Options:
	--seed=SEED
		Random seed for reproducibility

	--n_snps_rand=N_SNPS_RAND
		Number of SNPs modeled as a random effect

	--n_snps_fixed=N_SNPS_FIXED
		Number of SNPs modeled as a fixed effect

	--n_populations=N_POPULATIONS
		Number of sub-populations with distinct allele frequencies

	--rand_snp_freq=RAND_SNP_FREQ
		Frequency for SNPs modeled as random effects. If
		comma-separated list, frequency in each sub-population.

	--fixed_snp_freq=FIXED_SNP_FREQ
		Frequency for SNPs modeled as fixed effects. If comma-separated
		list, frequency in each sub-population.

	--n_samples=N_SAMPLES
		Number of samples. If comma-separated list, number of samples
		in each sub-population.

	--beta_fixed=BETA_FIXED
		Comma-separated list of beta values for the fixed effect SNPs.
		Must be of length equal to 'n_snps_fixed'

	--intercept=INTERCEPT
		Intercept of the phenotype

	--ploidy=PLOIDY
		Ploidy of the organism. Affects the range of dosages, which
		will be [0,ploidy]

	--heritability=HERITABILITY
		Narrow-sense heritability of the simulated phenotype. Must be
		between 0 and 1.

	--variance_scale=VARIANCE_SCALE
		Scale of the genetic and environmental variances. The genetic
		variance is heritability * variance_scale, while the
		environmental variance is (1 - heritability) * variance_scale

	--n_qcov=N_QCOV
		Number of quantitative covariates

	--beta_qcov=BETA_QCOV
		Beta for the quantitative covariate. If comma-separated list,
		beta for each of the quantitative covariates. Must be of length
		'n_qcov'

	--mean_qcov=MEAN_QCOV
		Mean for the quantitative covariate. If comma-separated list,
		mean for each of the quantitative covariates. Must be of length
		'n_qcov'

	--sd_qcov=SD_QCOV
		Standard deviation for the quantitative covariate. If
		comma-separated list, standard deviation for each of the
		quantitative covariates. Must be of length 'n_qcov'

	--n_cov=N_COV
		Number of covariates

	--n_levels_cov=N_LEVELS_COV
		Number of levels for the covariate. If comma-separated list,
		number of levels for each of the covariates. Must be of length
		'n_cov'. Each covariate must have > 1 level

	--beta_cov=BETA_COV1_LEV1,BETA_COV1_LEV2,BETA_COV2_LEV1,BETA_COV2_LEV2
		Beta for each level of the covariate. Comma-separated list with
		one value for each level of the covariate - 1 (dummy encoding
		is used), with betas for the levels of multiple covariates
		given sequentially. Must be of length equal to the sum of
		(n_levels - 1) for each covariate.

	--out=OUT
		Output name prefix

	-h, --help
		Show this help message and exit
```
