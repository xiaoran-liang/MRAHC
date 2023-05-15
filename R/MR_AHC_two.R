

#' Implementing the MR-AHC method for detecting genetic variants clusters
#' and conducting post-clustering estimation using summary statistics in Mendelian randomization with two outcomes.

#'@param betaX: A numeric vector of SNP-exposure associations (p by 1 vector).
#'@param betaY1: A numeric vector of SNP-outcome associations for outcome 1 (p by 1 vector).
#'@param betaY2: A numeric vector of SNP-outcome associations for outcome 2 (p by 1 vector).
#'@param seX: A numeric vector of the standard errors of the SNP-exposure associations (p by 1 vector).
#'@param seY1: A numeric vector of the standard errors of the SNP-outcome associations for outcome 1 (p by 1 vector).
#'@param seY2: A numeric vector of the standard errors of the SNP-outcome associations for outcome 2 (p by 1 vector).
#'@param n: A numeric scalar specifying the sample size. In a multi-sample MR design, we recommend using the smaller sample size of the outcome samples.
#'@param alpha: A numeric scalar between 0 and 1 specifying the significance level for the confidence interval (default = 0.05).
#'@param tuning: A numeric scalar specifying the threshold p-value for the Q test (default = 0.1/log(n)).
#'@param smallcluster: A numeric scalar specifying the cutoff value for small clusters (default = 4).
#'@param outremove: Logical. If outremove = TRUE, conduct the outlier removal procedure (default = FALSE).
#'@param iter: Logical. When outremove = TRUE, if iter = TRUE, conduct the iterated outlier removal procedure (default = FALSE).
#'@param iter.p: A numeric scalar specifying the threshold p-value for the individual Q in the outlier removal procedure.
#'@param rho: A numeric scalar specifying the covariance between the two SNP-specific Wald estimates for the two outcomes.
#'            If the two outcome samples are non-overlapping, set rho = 0 (default = 0).

#'@return Cluster_number: The number of clusters detected by the algorithm (including the junk cluster).
#'@return Cluster_number_real: The number of substantive clusters detected by the algorithm (excluding the junk cluster).
#'@return AHC_cluster: All the clusters detected by the algorithm. SNPs identities (as in the dataset) are listed in each cluster.
#'@return AHC_cluster_real: All the substantive clusters detected by the algorithm. SNPs identities (as in the dataset) are listed in each cluster.
#'@return Null_cluster_1: Cluster index of the null cluster for outcome 1.
#'@return Null_cluster_2: Cluster index of the null cluster for outcome 2.
#'@return Junk_cluster: Identities of SNPs in the junk cluster.
#'@return F: the F statistic for all the IVs.
#'@return AHC_results: A matrix that summarizes the clustering and estimation results, including:
#'           (a) length: The number of IVs in each cluster.
#'           (b) beta: the point estimates estimated with each cluster.
#'           (c) se: the standard errors for the causal estimates in each cluster
#'           (d) t: the t statistic.
#'           (e) Qp: The p-value for the Q test of the instruments in each cluster.
#'           (f) I^2: The I statistic of the instruments in each cluster.
#'@return Confidence interval: The 95% (default) confidence intervals of the cluster-specific estimates.
#'@export

MR_AHC_two <-function(betaX, betaY1, betaY2, seX, seY1, seY2, n, alpha = 0.05,
                      tuning = 0.1/log(n), smallcluster = 4, outremove = FALSE,
                      iter = FALSE, iter.p = 0.05, rho = 0){

  # Check Input
  if(is.data.frame(betaX)){betaX <- as.matrix(betaX)}
  if(!is.vector(betaX) && !is.matrix(betaX) | !is.numeric(betaX) | ncol(as.matrix(betaX))!=1)
    stop("betaX must be a numeric vector.");
  betaX = as.numeric(betaX);

  if(is.data.frame(betaY1)){betaY1 <- as.matrix(betaY1)}
  if(!is.vector(betaY1) && !is.matrix(betaY1) | !is.numeric(betaY1) | ncol(as.matrix(betaY1))!=1)
    stop("betaY1 must be a numeric vector.");
  betaY1 = as.numeric(betaY1);

  if(is.data.frame(betaY2)){betaY2 <- as.matrix(betaY2)}
  if(!is.vector(betaY2) && !is.matrix(betaY2) | !is.numeric(betaY2) | ncol(as.matrix(betaY2))!=1)
    stop("betaY2 must be a numeric vector.");
  betaY2 = as.numeric(betaY2);

  if(is.data.frame(seX)){seX <- as.matrix(seX)}
  if(!is.vector(seX) && !is.matrix(seX) | !is.numeric(seX) | ncol(as.matrix(seX))!=1)
    stop("seX must be a numeric vector.");
  seX = as.numeric(seX);

  if(is.data.frame(seY1)){seY1 <- as.matrix(seY1)}
  if(!is.vector(seY1) && !is.matrix(seY1) | !is.numeric(seY1) | ncol(as.matrix(seY1))!=1)
    stop("seY1 must be a numeric vector.");
  seY1 = as.numeric(seY1);

  if(is.data.frame(seY2)){seY2 <- as.matrix(seY2)}
  if(!is.vector(seY2) && !is.matrix(seY2) | !is.numeric(seY2) | ncol(as.matrix(seY2))!=1)
    stop("seY2 must be a numeric vector.");
  seY2 = as.numeric(seY2);

  # Define Constants
  pz <- length(betaX); ## the number of SNPs

  # F statistic
  F_bar = mean(betaX^2/seX^2);

  # Wald ratio
  bvec1 <- matrix(betaY1/betaX, ncol = 1);
  bvec2 <- matrix(betaY2/betaX, ncol = 1);
  sevec1 <- matrix(abs(seY1/betaX), ncol = 1);
  sevec2 <- matrix(abs(seY2/betaX), ncol = 1);

  # Q test involving all instruments
  Qsq.all <- Q_two(sk = 1:pz, b1 = bvec1, b2 = bvec2, se1 = sevec1, se2 = sevec2, rho = rho);
  Qsq.all.p <- Qsq.all$Q_p;

  if (Qsq.all.p > tuning){
    # all the IVs are selected as in the same cluster. No need for clustering.
    AHC_real <- list(c(1:pz));
    AHC_junk <- NULL;
  }else{
    # AHC procedure
    AHC_cluster <- AHC.IV.two(bvec1, bvec2, sevec1, sevec2, tuning = tuning, rho = rho) ## this is the substantive clusters
    # IVs that do not belong to any clusters are classified as in the junk cluster.
    AHC_junk <- sort(c(1:pz)[-unlist(AHC_cluster)]);

    # trimming the small clusters
    small.index = which(lapply(AHC_cluster, length) <= smallcluster)
    if (length(small.index) != 0){
      AHC_junk_new <- sort(unlist(AHC_cluster[small.index]))
      AHC_junk = sort(union(AHC_junk_new, AHC_junk))
      AHC_cluster <- AHC_cluster[-small.index]
    }

    AHC_real <- AHC_cluster; #substantive clusters
    Nr_real <- length(AHC_real); # number of substantive clusters

    ## Outlier removal
    if (outremove == TRUE){

      ## Iterated outlier removal
      if (iter == TRUE){

        n_real <- c();
        n_real[1] <- Inf;
        n_real[2] <- length(unlist(AHC_real)); ## number of SNPs in the substantive clusters

        i = 2;
        while (n_real[i] - n_real[i - 1] < 0){

          Nr_real <- length(AHC_real);
          # filter SNPs with large contribution to the Q statistic
          for (k_iter in 1:Nr_real){
            leaf <- AHC_real[[k_iter]];
            Q_ind <- Q_two(sk = leaf, b1 = bvec1, b2 = bvec2, se1 = sevec1, se2 = sevec2, rho = rho)$Q_ind;
            Q_ind_p <- 1-pchisq(Q_ind, df = 2);

            p_threshold <- iter.p;
            leaf_new <- leaf[which(Q_ind_p >= p_threshold)];
            leaf_remove <- leaf[which(Q_ind_p < p_threshold)];

            AHC_real[[k_iter]] <- leaf_new;
            AHC_junk <- sort(union(AHC_junk, leaf_remove))
          }

          real_snp <- sort(unlist(AHC_real));
          i <- i + 1;
          n_real[i] <- length(real_snp);

          ## if any SNPs get removed, re-do clustering
          if (n_real[i] - n_real[i - 1] < 0){

            AHC_new <- AHC.IV.two(bvec1[real_snp], bvec2[real_snp], sevec1[real_snp], sevec2[real_snp], tuning = tuning, rho = rho);

            AHC_real <- AHC_new ## re-value AHC_real;
            for (k_new in 1:length(AHC_new)){
              AHC_real[[k_new]] <- real_snp[AHC_new[[k_new]]];
            }

            junk_new <- setdiff(real_snp, sort(unlist(AHC_real)));
            AHC_junk <- sort(c(AHC_junk, junk_new))

            # trimming the small clusters
            small.index = which(lapply(AHC_real, length) <= smallcluster)
            if (length(small.index) != 0){
              AHC_junk_new <- sort(unlist(AHC_real[small.index]))
              AHC_junk = sort(union(AHC_junk_new, AHC_junk))
              AHC_real <- AHC_real[-small.index]
            }
          }

        }## end of while
      }else{

        ## one-time outlier removal without iteration
        # filter SNPs with large contribution to the Q statistic
        n_real_init <- length(unlist(AHC_real));
        for (k_iter in 1:Nr_real){
          leaf <- AHC_real[[k_iter]];
          Q_ind <- Q_two(sk = leaf, b1 = bvec1, b2 = bvec2, se1 = sevec1, se2 = sevec2, rho = rho)$Q_ind;
          Q_ind_p <- 1-pchisq(Q_ind, df = 2);

          p_threshold <- iter.p;
          leaf_new <- leaf[which(Q_ind_p >= p_threshold)];
          leaf_remove <- leaf[which(Q_ind_p < p_threshold)];

          AHC_real[[k_iter]] <- leaf_new;
          AHC_junk <- sort(union(AHC_junk, leaf_remove))
        }
        real_snp <- sort(unlist(AHC_real));
        n_real <- length(real_snp);
        ## if any SNPs get removed, re-do clustering
        if (n_real < n_real_init){

          AHC_new <- AHC.IV.two(bvec1[real_snp], bvec2[real_snp], sevec1[real_snp], sevec2[real_snp], tuning = tuning, rho = rho);

          AHC_real <- AHC_new ## re-value AHC_real;
          for (k_new in 1:length(AHC_new)){
            AHC_real[[k_new]] <- real_snp[AHC_new[[k_new]]];
          }

          junk_new <- setdiff(real_snp, sort(unlist(AHC_real)));
          AHC_junk <- sort(c(AHC_junk, junk_new))

          # trimming the small clusters
          small.index = which(lapply(AHC_real, length) <= smallcluster)
          if (length(small.index) != 0){
            AHC_junk_new <- sort(unlist(AHC_real[small.index]))
            AHC_junk = sort(union(AHC_junk_new, AHC_junk))
            AHC_real <- AHC_real[-small.index]
          }
        }

      }

    } ## end of outlier removal
  }


  AHC_cluster = AHC_cluster_real <- AHC_real;
  AHC_cluster_junk <- AHC_junk;

  if (length(AHC_cluster_junk) != 0){AHC_cluster[[length(AHC_cluster) + 1]] <- AHC_cluster_junk}else{AHC_cluster_junk <- NULL}

  Nr_real <- length(AHC_cluster_real); # number of substantive clusters
  Nr_all <- length(AHC_cluster) # number of total clusters including the junk cluster.

  # Post-clustering estimation
  AHC_results <- matrix(NA, nrow = Nr_real, ncol = 10);
  colnames(AHC_results) <- c("ID", "length", "beta1", "se1", "t1", "beta2", "se2", "t2", "Qp", "I^2");

  for (k_r in 1:Nr_real){
    AHC_results[k_r, "ID"] <- k_r;
    AHC_results[k_r, "length"] <- length(AHC_cluster_real[[k_r]]);

    leaf_kr <- AHC_cluster_real[[k_r]];
    results_kr <- Q_two(sk = leaf_kr, b1 = bvec1, b2 = bvec2, se1 = sevec1, se2 = sevec2, rho = rho)
    AHC_results[k_r, c("beta1", "beta2")] <- c(results_kr$ivw_est_k1, results_kr$ivw_est_k2);
    AHC_results[k_r, c("se1", "se2")] <- sqrt(c(results_kr$ivw_var_k1, results_kr$ivw_var_k2));
    AHC_results[k_r, c("t1", "t2")] <- abs(AHC_results[k_r, c("beta1", "beta2")]/AHC_results[k_r, c("se1", "se2")]);
    AHC_results[k_r, "Qp"] <- results_kr$Q_p;
    AHC_results[k_r, "I^2"] <- results_kr$Isq
  }

  confidence_interval <- matrix(NA, nrow = Nr_real, ncol = 4);
  rownames(confidence_interval) <- paste0("Cluster", 1:Nr_real);
  colnames(confidence_interval) <- c("lower-beta1", "upper-beta1", "lower-beta2", "upper-beta2");
  confidence_interval[, c("lower-beta1", "lower-beta2")] <-  AHC_results[, c("beta1", "beta2")]-abs(qnorm(alpha/2))*AHC_results[, c("se1", "se2")];
  confidence_interval[, c("upper-beta1", "upper-beta2")] <-  AHC_results[, c("beta1", "beta2")]+abs(qnorm(alpha/2))*AHC_results[, c("se1", "se2")];

  # The null clusters
  null_index_1 <- which(pchisq((AHC_results[,"t1"])^2, df = 1, lower.tail = FALSE) > tuning);
  null_index_2 <- which(pchisq((AHC_results[,"t2"])^2, df = 1, lower.tail = FALSE) > tuning);

  # Report results
  results = list( Cluster_number = Nr_all,
                  Cluster_number_real = Nr_real,
                  AHC_cluster = AHC_cluster,
                  AHC_cluster_real = AHC_cluster_real,
                  Null_cluster_1 = null_index_1,
                  Null_cluster_2 = null_index_2,
                  Junk_cluster = AHC_cluster_junk,
                  F = F_bar,
                  AHC_results = AHC_results,
                  confidence_interval = confidence_interval)
  return(results)

}


## below are internal functions, not to be called by users.

## IVW
ivw_est <- function(b, se){
  w <- 1/se^2;
  est_ivw <- sum(b * w)/sum(w);
  var_ivw <- 1/sum(w);
  return(c(est_ivw, var_ivw))
}

## Wald test
wald_test <- function(b1, se1, b2, se2){
  return((b1 - b2)^2 / (se1^2 + se2^2))
}

## D_kl
dissimi <- function(sk, sl, b1, b2, se1, se2, rho){

  bk1 <- b1[sk];
  sek1 <- se1[sk];
  bk2 <- b2[sk];
  sek2 <- se2[sk];

  bl1 <- b1[sl];
  sel1 <- se1[sl];
  bl2 <- b2[sl];
  sel2 <- se2[sl];

  ivw_k1 <- ivw_est(bk1, sek1);
  ivw_k2 <- ivw_est(bk2, sek2);
  ivw_l1 <- ivw_est(bl1, sel1);
  ivw_l2 <- ivw_est(bl2, sel2);

  ivw_est_k1 <- ivw_k1[1];
  ivw_var_k1 <- ivw_k1[2];
  ivw_est_k2 <- ivw_k2[1];
  ivw_var_k2 <- ivw_k2[2];

  ivw_est_l1 <- ivw_l1[1];
  ivw_var_l1 <- ivw_l1[2];
  ivw_est_l2 <- ivw_l2[1];
  ivw_var_l2 <- ivw_l2[2];

  ivw_diff <- matrix(c(ivw_est_k1 - ivw_est_l1, ivw_est_k2 - ivw_est_l2), ncol = 1);

  var_kl1 <- ivw_var_k1 + ivw_var_l1;
  var_kl2 <- ivw_var_k2 + ivw_var_l2;
  cov_kl <- rho * ivw_var_k1 * ivw_var_k2 * sum(1/sek1^2 * 1/sek2^2) +
    rho * ivw_var_l1 * ivw_var_l2 * sum(1/sel1^2 * 1/sel2^2);

  var_matrix <- matrix(c(var_kl1, cov_kl, cov_kl, var_kl2), nrow = 2);

  d_kl <- t(ivw_diff) %*% solve(var_matrix) %*% ivw_diff;

  return(list(dkl = d_kl,
              ivw_ests = matrix(c(ivw_est_k1, ivw_est_k2, ivw_est_l1, ivw_est_l2), ncol = 2),
              var_matrix = var_matrix))

}

## merging process
hier_tree <- function(b1, b2, se1, se2, rho){

  tree_id <- list();
  tree_est <- list();
  tree_dissm <- list();
  tree_len <- length(b1);

  tree_id[[1]] <- as.list(c(1:tree_len));
  tree_est[[1]] <- cbind(b1, b2, se1, se2);

  ## initial step
  base_id_init <- tree_id[[1]];
  base_est_init <- tree_est[[1]];
  dissm_init <- matrix(NA, nrow = tree_len, ncol = tree_len);

  for (j in 1:(tree_len - 1)){
    for (k in (j + 1):tree_len){
      dissm_init[j, k] <- dissimi(base_id_init[[j]], base_id_init[[k]], b1, b2, se1, se2, rho = 0)[[1]];
    }
  }

  tree_dissm[[1]] <- dissm_init;

  merge_id_init <- c(which(dissm_init == min(dissm_init, na.rm = TRUE), arr.ind = TRUE));
  merge_est_init <- c(ivw_est(b1[merge_id_init], se1[merge_id_init])[1],
                      ivw_est(b2[merge_id_init], se2[merge_id_init])[1],
                      sqrt(ivw_est(b1[merge_id_init], se1[merge_id_init])[2]),
                      sqrt(ivw_est(b2[merge_id_init], se2[merge_id_init])[2]));

  base_id <- base_id_init[-merge_id_init];
  base_id[[length(base_id) + 1]] <- merge_id_init;

  base_est <- base_est_init[-merge_id_init, ];
  base_est <- rbind(base_est, merge_est_init);

  base_dissm <- dissm_init[-merge_id_init, -merge_id_init];
  #base_dissm <- cbind(base_dissm, rep(NA, nrow(base_dissm))) %>% rbind(rep(NA, ncol(base_dissm)+1));
  base_dissm <- rbind(cbind(base_dissm, rep(NA, nrow(base_dissm))),
                      rep(NA, ncol(base_dissm)+1));

  tree_id[[2]] <- base_id;
  tree_est[[2]] <- base_est;
  tree_dissm[[2]] <- base_dissm;

  i = 2;
  while (i < tree_len - 1) {

    base_id <- tree_id[[i]];
    base_est <- tree_est[[i]];
    base_dissm <- tree_dissm[[i]];

    branch_length <- length(base_id);

    ## only need to calculate the dissimilarity for the new cluster
    diss_temp <- rep(NA, branch_length - 1);
    for (j in 1:(branch_length - 1)){
      diss_temp[j] <- dissimi(base_id[[branch_length]], base_id[[j]], b1, b2, se1, se2, rho = 0)[[1]];
      tree_dissm[[i]][j, branch_length] <- diss_temp[j];
    }

    merge_id <- c(which(tree_dissm[[i]] == min(tree_dissm[[i]], na.rm = TRUE), arr.ind = TRUE));

    #if (length(merge_id) > 2){merge_id <- merge_id[c(1,3)]}

    merge_id_ind <- sort(unlist(base_id[merge_id]));

    merge_est <- c(ivw_est(b1[merge_id_ind], se1[merge_id_ind])[1],
                   ivw_est(b2[merge_id_ind], se2[merge_id_ind])[1],
                   sqrt(ivw_est(b1[merge_id_ind], se1[merge_id_ind])[2]),
                   sqrt(ivw_est(b2[merge_id_ind], se2[merge_id_ind])[2]));

    i <- i + 1;

    tree_id[[i]] <- base_id[-merge_id];
    tree_id[[i]][[branch_length - 1]] <- merge_id_ind;

    tree_est[[i]] <- base_est[-merge_id, ];
    tree_est[[i]] <- rbind(tree_est[[i]], merge_est);

    tree_dissm[[i]] <- as.matrix(tree_dissm[[i - 1]][-merge_id, -merge_id]);
    #tree_dissm[[i]] <- cbind(tree_dissm[[i]], rep(NA, nrow(tree_dissm[[i]]))) %>% rbind(rep(NA, ncol(tree_dissm[[i]])+1));
    tree_dissm[[i]] <- rbind(cbind(tree_dissm[[i]], rep(NA, nrow(tree_dissm[[i]]))),
                             rep(NA, ncol(tree_dissm[[i]])+1));

  }

  merge_id_ind <- c(1:tree_len);
  tree_id[[tree_len]] <- merge_id_ind;
  tree_est[[tree_len]] <- c(ivw_est(b1[merge_id_ind], se1[merge_id_ind])[1],
                            ivw_est(b2[merge_id_ind], se2[merge_id_ind])[1],
                            sqrt(ivw_est(b1[merge_id_ind], se1[merge_id_ind])[2]),
                            sqrt(ivw_est(b2[merge_id_ind], se2[merge_id_ind])[2]));

  return(list(tree_id = tree_id,
              tree_est = tree_est,
              tree_dissm = tree_dissm))

}


## two_outcome Q
Q_two <- function(sk, b1, b2, se1, se2, rho){

  sk_len <- length(sk);
  df <- 2*(sk_len - 1);

  bk1 <- b1[sk];
  sek1 <- se1[sk];
  bk2 <- b2[sk];
  sek2 <- se2[sk];

  ivw_k1 <- ivw_est(bk1, sek1);
  ivw_k2 <- ivw_est(bk2, sek2);

  ivw_est_k1 <- ivw_k1[1];
  ivw_est_k2 <- ivw_k2[1];

  ivw_var_k1 <- ivw_k1[2];
  ivw_var_k2 <- ivw_k2[2];

  var_matrix <- rbind(cbind(diag(sek1^2), diag(rho, sk_len)),
                      cbind(diag(rho, sk_len), diag(sek2^2)));

  Q_crit <- matrix(c(bk1 - rep(ivw_est_k1, sk_len), bk2 - rep(ivw_est_k2, sk_len)), nrow = 1) %*%
    solve(var_matrix) %*%
    matrix(c(bk1 - rep(ivw_est_k1, sk_len), bk2 - rep(ivw_est_k2, sk_len)), ncol = 1);

  Q_ind <- rep(NA, sk_len)
  for (k in 1:sk_len){

    var_matrix_ind <- rbind(cbind(diag(sek1[k]^2, 1), diag(rho, 1)),
                            cbind(diag(rho, 1), diag(sek2[k]^2, 1)));

    Q_ind_value <- matrix(c(bk1[k] - rep(ivw_est_k1, 1), bk2[k] - rep(ivw_est_k2, 1)), nrow = 1) %*%
      solve(var_matrix_ind) %*%
      matrix(c(bk1[k] - rep(ivw_est_k1, 1), bk2[k] - rep(ivw_est_k2, 1)), ncol = 1);

    #Q_ind[k] <- Q_ind_value/Q_crit
    Q_ind[k] <- Q_ind_value
  }

  Q_p <- 1-pchisq(Q_crit,df);

  Isq <- (Q_crit - df)/Q_crit
  Isq <- max(0,Isq)

  return(list(Q_crit = Q_crit, Q_p = Q_p, ivw_est_k1 = ivw_est_k1, ivw_est_k2 = ivw_est_k2,
              ivw_var_k1 = ivw_var_k1, ivw_var_k2 = ivw_var_k2, Isq = Isq, Q_ind = Q_ind))
}

# Internal Function, not to be called by users.
AHC.IV.two <- function(bvec1, bvec2, sevec1, sevec2, tuning, rho){

  # Define Constants
  pz <- length(bvec1); ## the number of SNPs

  # AHC procedure
  hiertree <- hier_tree(bvec1, bvec2, sevec1, sevec2, rho)$tree_id;

  # Downward testing
  start <- length(hiertree) - 1;
  clusterIV = Q.IV = beta.IV = se.IV = I.IV = clusternames <- list();
  while (start > 1) {
    branch = hiertree[[start]];
    index_single = which(lapply(branch, length) == 1);
    if (length(index_single) != 0) branch = branch[-index_single];
    for (k in 1:length(branch)){
      leaf <- sort(unique(as.vector(branch[[k]])));
      check.subset <- function(y){all(leaf %in% y)};
      index = any(unlist(lapply(clusterIV, check.subset)));
      if (length(leaf)< pz && index == FALSE){
        Q.branch <- Q_two(leaf, bvec1, bvec2, sevec1, sevec2, rho)[[1]];
        p.branch <- Q_two(leaf, bvec1, bvec2, sevec1, sevec2, rho)[[2]];
        beta.branch <- unlist(Q_two(leaf, bvec1, bvec2, sevec1, sevec2, rho)[3:4]);
        se.branch <- sqrt(unlist(Q_two(leaf, bvec1, bvec2, sevec1, sevec2, rho)[5:6]));
        I.branch <- Q_two(leaf, bvec1, bvec2, sevec1, sevec2, rho)[[7]];
        if (p.branch > tuning){
          #if (I.branch < 0.2){
          clusterIV[[length(clusterIV)+1]] <- leaf;
          Q.IV[[length(Q.IV)+1]] <- p.branch;
          beta.IV[[length(beta.IV)+1]] <- beta.branch;
          se.IV[[length(se.IV)+1]] <- se.branch;
          I.IV[[length(se.IV)+1]] <- I.branch;
          clusternames[[length(clusternames)+1]] <- paste0('cluster', length(clusternames)+1)
        }
      }
    }
    clusterIV;
    start <- start - 1;
  }

  return(clusterIV)
}
