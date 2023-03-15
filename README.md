# MR-AHC
MR-AHC is an R package for detecting genetic variant clusters and conducting post-clustering estimation using summary statistics in Mendelian randomisation. When there is mechanistic homogeneity, genetic variants involved in MR analysis can be divided into clusters identifying distinct causal effects driven by different biological mechanisms. Variants are classified into the same cluster only if their variant-specific causal estimates are similar to each other.

# Functions
There are two main functions:

- **MR_AHC**: conduct clustering analysis in MR with one outcome and one exposure.
- **MR_AHC_two**: conduct clustering analysis in MR with two outcomes and a common exposure.

# Installation

```R
install.packages("devtools")
devtools::install_github("xiaoran-liang/MR-AHC")
```

# Running example
```R
## The usage of the two main functions are illustrated with simulated data generated with sample size = 60000
## (1). illustration of the usage of the MR_AHC function
## read data
dat <- MRAHC::date_example_one
betaX <- dat$est_x  # SNP-exposure associations
seX <- dat$se_x  # standard errors of the SNP-exposure associations
betaY <- dat$est_y  # SNP-outcome associations
seY <- dat$se_y  # standard errors of the SNP-outcome associations

## run the function without outlier removal
output_1_a <- MR_AHC(betaX, betaY, seX, seY, n = 60000, outremove = FALSE)

## run the function with one-time outlier removal
output_1_b <- MR_AHC(betaX, betaY, seX, seY, n = 60000, outremove = TRUE, iter = FALSE)

## run the function with iterated outlier removal
output_1_c <- MR_AHC(betaX, betaY, seX, seY, n = 60000, outremove = TRUE, iter = TRUE)

## check the output. All three ways have the same output, take the one without outlier removal as an example:
## 1. total number of detected clusters
output_1_a$Cluster_number
## 2. the number of non-junk clusters
output_1_a$Cluster_number_real
## 3. genetic variant clusters detected by the algorithm
output_1_a$AHC_cluster
## 4. non-junk genetic variant clusters detected by the algorithm
output_1_a$AHC_cluster_real
## 5. the cluster index of the null cluster
output_1_a$Null_cluster
## 6. genetic variants in the junk cluster
output_1_a$Junk_cluster
## 7. the F statstic
output_1_a$F
## 8. the post-clustering estimation results, including:
##           (a) length: The number of IVs in each cluster.
##           (b) beta: the point estimate estimated with each cluster.
##           (c) se: the standard error for the causal estimate in each cluster
##           (d) t: the t statistic.
##           (e) Qp: The p-value for the Q test of the instruments in each cluster.
##           (f) I^2: The I statistic of the instruments in each cluster.
output_1_a$AHC_results
## 9. the 95% confidence intervals for the cluster-specific IVW estimates
output_1_a$confidence_interval

## (2). illustration of the usage of the MR_AHC_two function
## read data
dat <- MRAHC::date_example_two
betaX <- dat$est_x  # SNP-exposure associations
seX <- dat$se_x  # standard errors of the SNP-exposure associations
betaY1 <- dat$est_y1  # SNP-outcome associations for outcome 1
seY1 <- dat$se_y1  # standard errors of the SNP-outcome associations for outcome 1
betaY2 <- dat$est_y2  # SNP-outcome associations for outcome 2
seY2 <- dat$se_y2  # standard errors of the SNP-outcome associations for outcome 2

## run the function without outlier removal
output_2_a <- MR_AHC_two(betaX, betaY1, betaY2, seX, seY1, seY2, n = 60000, outremove = FALSE)

## run the function with one-time outlier removal
output_2_b <- MR_AHC_two(betaX, betaY1, betaY2, seX, seY1, seY2, n = 60000, outremove = TRUE, iter = FALSE)

## run the function with iterated outlier removal
output_2_c <- MR_AHC_two(betaX, betaY1, betaY2, seX, seY1, seY2, n = 60000, outremove = TRUE, iter = TRUE)

## check the output. The output is the same as the MR_AHC function, with the exception of the null cluster results:
output_2_c$Null_cluster_1 ## the cluster index of the null cluster for outcome 1
output_2_c$Null_cluster_2 ## the cluster index of the null cluster for outcome 2

```

