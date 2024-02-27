########## make a set of snps that adhere to HWE
# init empty df
nobs <- 10000
nsnps <- 10
snps <- data.frame(matrix(0, nrow=nobs, ncol=nsnps))

# draw maf from uniform distribution
maf <- runif(nsnps, 0.35, 0.5)

# draw snp from multinomial distribution for each maf
for (i in 1:length(maf)) {
  snps[, i] <- rbinom(nobs, 2, maf[i])
}

######### derive the outcome
# set m_betas, the main effect sizes for SNPs
# all null for now, but could be e.g. rnorm(nobs, 0, 0.1)
# can check by calling exp() on the betas and seeing if they're reasonable odds ratios
# note that this is the effect of the SNPs on the liability scale (the underlying risk), not
# on the outcome scale (disease status), so ORs can be higher than a typically GWAS reports
m_betas <- rep(0, nsnps)

# set i_betas, the interaction effect sizes for SNPs
i_betas <- c(0.01, 0.01, 0, 0, 0)

# set which SNPs will interact
# currently only SNPs 1/2 and SNPs 3/4 will interact
# betas for these are above
i_combinations <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10))

# get g, the genotypic score
# this is the dot product of the genotypes and the effect sizes
# i.e. sum of each snp times is effect
# gives a score for each individual
g_main <- as.matrix(snps) %*% m_betas

# combine interaction effects under a multiplicative model
g_interaction <- rep(0, nobs)
for (i in 1: length(i_combinations)) {
  g_interaction <- g_interaction + ((snps[, i_combinations[[i]][1]] * snps[, i_combinations[[i]][2]]) * i_betas[i])
}

# sum together main and interaction effects
g <- g_main + g_interaction

# setting noise, or environmental risk, e
# here we assume all noise is environmental risk and lum into one risk score
# this is all the risk influecing underlying risk which we don't measure using genetics
# note that I've chosen something manually here, but in practice we should
# set this using the heritability or something else to make the signal/noise ratio realistic
e <- rnorm(nobs, 0, 0.01)

# our liability for disease, l, is the result of our genetic risk, g, plus our environmental risk, e
# where the genetic risk incorporates main and interaction effects
# note that here we just sum the two, so there are no snp-environment interactions
l <- g + e

# setting the threshold for disease status
# above this score we assign disease status to participants
# note that the threshold must be applied to the final liability, l, rather than the genotypic risk, g,
# otherwise we would get perfect prediction from genetics.
# we base this on the prevalence, k, of the disease in the population
# note that the value for prevalence here would make it a *very* common disease so is unrealistic for most scenarios
k <- 0.20

# we standardise the liability scale so we get the appropriate percentile
l_z <- (l - mean(l)) / sd(l)

# we then take that percentile of the standard normal distribution
threshold <- qnorm(1-k)

# phenotype is 1 (case) if risk exceeds threshold for disease, else 0 (control)
p <- as.integer(l_z > threshold)

# checking case/control ratio
# you may want to sub-sample controls to get 50-50 cases and controls,
# but it's easier to use reweighting when fitting classifiers
paste('Proportion of cases in sample: ', mean(p), sep='')

# combination into final df for checking
snps$p <- p

# checking with a glm and main effects only
m1 <- glm('p ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10',
         data=snps, family = 'binomial')
print(summary(m1))

# checking with interaction terms
m2 <- glm('p ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X1:X2 + X3:X4',
          data=snps, family = 'binomial')
print(summary(m2))

# with interaction terms after scaling
snps_z <- as.data.frame(scale(snps[,1:nsnps]))
snps_z$p <- snps$p

m3 <- glm('p ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X1:X2 + X3:X4',
          data=snps_z, family = 'binomial')
print(summary(m3))

