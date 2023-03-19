
# one outcome
launch_simulations <- function(N_x,    # exposure sample size
                               N_y,        # outcome sample size
                               SNPs_group, # matrix defining groups (different components / pleiotropy)
                               beta,       # causal effects of each sub-component of the exposure on outcome) (list of size K)
                               # direct effects of SNPs on X are generated from a uniform distribution
                               GX_u, # specify the upper limit of the uniform distribution
                               GX_l, # specify the lower limit of the uniform distribution
                               h2_K,      # heritability for the different component of exposure  (vector of size K)
                               # effects through a confounder (correlated pleiotropy)
                               # if no pleiotropy keep the following parameters = 0
                               h2_U = 0,  # heritability for the confounder
                               q_x = 0,   # effect of the confounder on the exposure
                               q_y =0   # effect of the confounder on the outcome
){
  
  # check if there's pleiotropy
  pleiotropy = F;
  if(!any(c(h2_U, q_x, q_y) == 0)) pleiotropy = T;
  
  K = nrow(SNPs_group); # number of groups
  M = ncol(SNPs_group); # number of SNPs
  MK = rowSums(SNPs_group)[1:K]; # number of SNPs in each group
  
  beta_vec <- rep(NA, M);
  check_y <- rep(NA, K);
  for (k in 1:K){
    SNP_k = which(SNPs_group[k, ] == 1)
    beta_vec[SNP_k] = beta[[k]];
    check_y[k] <- sum(h2_K[k]/MK[k] * beta_vec[SNP_k]**2);
  }
  check_y <- sum(check_y);
  
  if(sum(h2_K)+h2_U*q_x**2>=1) stop("X can not have a variance of 1 if sum(h2_K) + h2_U * q_x**2 >= 1.");
  if(check_y+h2_U*q_y**2>=1) stop("Y can not have a variance of 1 if  sum(h2_K * beta**2) + h2_U * q_y**2 >= 1.");
  
  # direct effects of SNPs on exposure
  effects_x = matrix(sample(c(-1,1), M, replace = T)*runif(M, GX_l, GX_u), 1, M);
  for(k in 1:K){
    SNPs_k = which(SNPs_group[k,]==1);
    effects_x[SNPs_k] = sqrt(h2_K[k])*effects_x[SNPs_k]/(sqrt(sum(effects_x[SNPs_k]^2))); # scale the SNP-exposure effects
  }
  
  # if there's pleiotropy, pleiotropic SNPs are set to be the first group
  M_U = sum(SNPs_group[1,]);
  SNP_U = which(SNPs_group[1,] == 1);
  effects_pleiotropy = matrix(0, 1, M);
  if(pleiotropy){
    effects_U = sample(c(-1,1), M_U, replace = T)*runif(M_U, 0.1, 0.3);
    effects_U = sqrt(h2_U)*effects_U/(sqrt(sum(effects_U^2)));
    effects_pleiotropy[SNP_U] = effects_U;
    effects_x[SNP_U] <- 0
  } else {
    pi_U = h2_U = q_x = q_y = 0;
  }
  
  # the SNP-exposure associations
  betaX = effects_x + effects_pleiotropy * q_x + rnorm(M, 0, sqrt(1/N_x));
  betaX_se = rep(1/sqrt(N_x), M);
  
  # the SNP-outcome associations
  betaY = effects_x * beta_vec + effects_pleiotropy * beta_vec * q_x  + effects_pleiotropy * q_y + rnorm(M, 0, sqrt(1/N_y));
  betaY_se = rep(1/sqrt(N_y), M);
  
  # the F statistic
  F_bar = mean(betaX^2/betaX_se^2);
  
  results = list( summary_statistics = data.frame( est_x = as.numeric(betaX), 
                                                   se_x = as.numeric(betaX_se), 
                                                   est_y = as.numeric(betaY), 
                                                   se_y = as.numeric(betaY_se)),
                  F = F_bar)
  return(results)
}


# two outcomes
launch_simulations_two <- function(N_x,    # exposure sample size
                                   N_y1,        # outcome 1 sample size
                                   N_y2,        # outcome 2 sample size
                                   SNPs_group, # matrix defining groups (different components / pleiotropy)
                                   beta1,       # causal effects of each sub-component of the exposure on outcome 1) (list of size K)
                                   beta2,       # causal effects of each sub-component of the exposure on outcome 2) (list of size K)
                                   # direct effects of SNPs on X are generated from a uniform distribution
                                   GX_u, # specify the upper limit of the uniform distribution
                                   GX_l, # specify the lower limit of the uniform distribution
                                   h2_K,      # heritability for the different component of exposure  (vector of size K)
                                   # effects through a confounder (correlated pleiotropy)
                                   # if no pleiotropy keep the following parameters = 0
                                   h2_U = 0,  # heritability for the confounder
                                   q_x = 0,   # effect of the confounder on the exposure
                                   q_y1 =0,   # effect of the confounder on the outcome 1
                                   q_y2 =0   # effect of the confounder on the outcome 2
){
  
  # check if there's pleiotropy
  pleiotropy = F;
  if(!any(c(h2_U, q_x, q_y1, q_y2) == 0)) pleiotropy = T;
  
  K = nrow(SNPs_group); # number of groups
  M = ncol(SNPs_group); # number of SNPs
  MK = rowSums(SNPs_group)[1:K]; # number of SNPs in each group
  
  beta_vec1 <- beta_vec2 <- rep(NA, M);
  check_y1 <- check_y2 <- rep(NA, K);
  for (k in 1:K){
    SNP_k = which(SNPs_group[k, ] == 1)
    beta_vec1[SNP_k] = beta1[[k]];
    beta_vec2[SNP_k] = beta2[[k]];
    check_y1[k] <- sum(h2_K[k]/MK[k] * beta_vec1[SNP_k]**2);
    check_y2[k] <- sum(h2_K[k]/MK[k] * beta_vec2[SNP_k]**2);
  }
  check_y1 <- sum(check_y1);
  check_y2 <- sum(check_y2);
  
  if(sum(h2_K)+h2_U*q_x**2>=1) stop("X can not have a variance of 1 if sum(h2_K) + h2_U * q_x**2 >= 1.");
  if(check_y1+h2_U*q_y1**2>=1) stop("Y can not have a variance of 1 if  sum(h2_K * beta**2) + h2_U * q_y**2 >= 1.");
  if(check_y2+h2_U*q_y2**2>=1) stop("Y can not have a variance of 1 if  sum(h2_K * beta**2) + h2_U * q_y**2 >= 1.");
  
  # direct effects of SNPs on exposure
  effects_x = matrix(sample(c(-1,1), M, replace = T)*runif(M, GX_l, GX_u), 1, M);
  for(k in 1:K){
    SNPs_k = which(SNPs_group[k,]==1);
    effects_x[SNPs_k] = sqrt(h2_K[k])*effects_x[SNPs_k]/(sqrt(sum(effects_x[SNPs_k]^2))); # scale the SNP-exposure effects
  }
  
  # if there's pleiotropy, pleiotropic SNPs are set to be the first group
  M_U = sum(SNPs_group[1,]);
  SNP_U = which(SNPs_group[1,] == 1);
  effects_pleiotropy = matrix(0, 1, M);
  if(pleiotropy){
    effects_U = sample(c(-1,1), M_U, replace = T)*runif(M_U, 0.1, 0.3);
    effects_U = sqrt(h2_U)*effects_U/(sqrt(sum(effects_U^2)));
    effects_pleiotropy[SNP_U] = effects_U;
    effects_x[SNP_U] <- 0
  } else {
    pi_U = h2_U = q_x = q_y1 = q_y2 = 0;
  }
  
  # the SNP-exposure associations
  betaX = effects_x + effects_pleiotropy * q_x + rnorm(M, 0, sqrt(1/N_x));
  betaX_se = rep(1/sqrt(N_x), M);
  
  # the SNP-outcome1 associations
  betaY1 = effects_x * beta_vec1 + effects_pleiotropy * beta_vec1 * q_x  + effects_pleiotropy * q_y1 + rnorm(M, 0, sqrt(1/N_y1));
  betaY1_se = rep(1/sqrt(N_y1), M);
  
  # the SNP-outcome2 associations
  betaY2 = effects_x * beta_vec2 + effects_pleiotropy * beta_vec2 * q_x  + effects_pleiotropy * q_y2 + rnorm(M, 0, sqrt(1/N_y2));
  betaY2_se = rep(1/sqrt(N_y2), M);
  
  # the F statistic
  F_bar = mean(betaX^2/betaX_se^2);
  
  results = list( summary_statistics = data.frame( est_x = as.numeric(betaX), 
                                                   se_x = as.numeric(betaX_se), 
                                                   est_y1 = as.numeric(betaY1), 
                                                   se_y1 = as.numeric(betaY1_se),
                                                   est_y2 = as.numeric(betaY2), 
                                                   se_y2 = as.numeric(betaY2_se)),
                  F = F_bar)
  return(results)
}


