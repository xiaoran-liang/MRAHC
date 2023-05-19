
#' Implementing the MR-AHC method for detecting genetic variants clusters
#' and conducting post-clustering estimation using summary statistics in Mendelian randomization with one outcome.

#'@param betaX: A numeric vector of the SNP-exposure associations (p by 1 vector).
#'@param betaY: A numeric vector of the SNP-outcome associations (p by 1 vector).
#'@param seX: A numeric vector of the standard errors of the SNP-exposure associations (p by 1 vector).
#'@param seY: A numeric vector of the standard errors of the SNP-outcome associations (p by 1 vector).
#'@param n: A numeric scalar specifying the sample size.
#'In a two-sample MR design, we recommend using the sample size of the outcome sample.
#'@param alpha: A numeric scalar between 0 and 1 specifying the significance level
#'for the confidence interval (default = 0.05).
#'@param tuning: A numeric scalar specifying the threshold p-value for the Q test (default = 0.1/log(n)).
#'@param weak: Logical. If weak = TRUE, the weak IV robust Q test is used (default = TRUE).
#'@param smallcluster: A numeric scalar specifying the cutoff value for small clusters (default = 4).
#'@param outremove: Logical. If outremove = TRUE, conduct the outlier removal procedure (default = FALSE).
#'@param iter: Logical. When outremove = TRUE, if iter = TRUE, conduct the iterated
#'outlier removal procedure (default = FALSE).
#'@param iter.p: A numeric scalar specifying the threshold p-value for the individual Q in
#' the outlier removal procedure.

#'@return Cluster_number: The number of clusters detected by the algorithm (including the junk cluster).
#'@return Cluster_number_real: The number of substantive clusters detected by the algorithm (excluding the junk cluster).
#'@return AHC_cluster: The cluster ids of all the SNPs. id = 0 means the SNPs are assigned to the junk cluster.
#'@return Null_cluster: Cluster index of the null cluster.
#'@return Junk_cluster: Identities of SNPs in the junk cluster.
#'@return F: the F statistic for all the IVs.
#'@return AHC_results: A matrix that summarizes the clustering and estimation results, including:
#'           (a) length: The number of IVs in each cluster.
#'           (b) beta: the point estimate estimated with each cluster.
#'           (c) se: the standard error for the causal estimate in each cluster
#'           (d) t: the t statistic.
#'           (e) Qp: The p-value for the Q test of the instruments in each cluster.
#'           (f) I^2: The I statistic of the instruments in each cluster.
#'@return Confidence interval: The 95% (default) confidence intervals of the cluster-specific estimates.
#'@export

MR_AHC <-function(betaX, betaY, seX, seY, n, alpha = 0.05,
                  tuning = 0.1/log(n), weak = TRUE, smallcluster = 4, outremove = FALSE,
                  iter = FALSE, iter.p = 0.05){

  # Check Input
  if(is.data.frame(betaX)){betaX <- as.matrix(betaX)}
  if(!is.vector(betaX) && !is.matrix(betaX) | !is.numeric(betaX) | ncol(as.matrix(betaX))!=1)
    stop("betaX must be a numeric vector.");
  betaX = as.numeric(betaX);

  if(is.data.frame(betaY)){betaY <- as.matrix(betaY)}
  if(!is.vector(betaY) && !is.matrix(betaY) | !is.numeric(betaY) | ncol(as.matrix(betaY))!=1)
    stop("betaY must be a numeric vector.");
  betaY = as.numeric(betaY);

  if(is.data.frame(seX)){seX <- as.matrix(seX)}
  if(!is.vector(seX) && !is.matrix(seX) | !is.numeric(seX) | ncol(as.matrix(seX))!=1)
    stop("seX must be a numeric vector.");
  seX = as.numeric(seX);

  if(is.data.frame(seY)){seY <- as.matrix(seY)}
  if(!is.vector(seY) && !is.matrix(seY) | !is.numeric(seY) | ncol(as.matrix(seY))!=1)
    stop("seY must be a numeric vector.");
  seY = as.numeric(seY);

  # Define Constants
  pz <- length(betaY);
  P <- 1; ## this is for single outcome--single exposure case

  # F statistic
  F_bar = mean(betaX^2/seX^2);

  # Functions for the Q Test
  if (weak){
    # weak-IV robust Q
    Qsq = function(bx,by,sx,sy){
      biv        = by/bx
      k          = length(bx)
      df         = k - 1
      PL2 = function(a){
        b = a[1]
        w = 1/(sy^2/bx^2 + (b^2)*sx^2/bx^2)
        q = sum(w*(biv - b)^2)
      }
      Bhat = optimize(PL2,interval=c(-10,10))$minimum
      W = 1/(sy^2/bx^2 + (Bhat^2)*sx^2/bx^2)
      QE = sum(W*(biv - Bhat)^2)
      QE_ind = W*(biv - Bhat)^2
      QEp = 1-pchisq(QE,df)
      Isq = (QE - df)/QE
      Isq =  max(0,Isq)
      sehat = sqrt(solve(sum(W)))
      return(list(QE, Bhat, sehat, Isq, QE_ind))
    }
  }else{
    Qsq = function(bx,by,sx,sy){
      k          = length(bx)
      df         = k - 1
      y          = matrix(by/bx, nrow = k)
      s          = matrix(sy/bx, nrow = k)
      w          = 1/s^2; sum.w  = sum(w)
      mu.hat     = sum(y*w)/sum.w
      Q          = sum(w*(y-mu.hat)^2)
      se.hat     = sqrt(solve(sum.w))
      Isq        = (Q - df)/Q
      Isq        =  max(0,Isq)
      return(list(Q, mu.hat, se.hat, Isq))
    }
  }

  # Q test involving all instruments
  Qsq.all <- Qsq(betaX, betaY, seX, seY);
  Qsq.all.p <- pchisq(Qsq.all[[1]], (pz - P), lower.tail = FALSE);

  if (Qsq.all.p > tuning){
    # all the IVs are selected as in the same cluster. No need for clustering.
    AHC_real <- list(c(1:pz));
    AHC_junk <- NULL;
  }else{
    # AHC procedure
    AHC_real <- AHC.IV(betaX, betaY, seX, seY, n, tuning, weak); ## this is the substantive clusters
    # IVs that do not belong to any clusters are classified as in the junk cluster.
    AHC_junk <- sort(c(1:pz)[-unlist(AHC_real)]);

    # trimming the small clusters
    small.index = which(lapply(AHC_real, length) <= smallcluster)
    if (length(small.index) != 0){
      AHC_junk_new <- sort(unlist(AHC_real[small.index]))
      AHC_junk = sort(union(AHC_junk_new, AHC_junk))
      AHC_real <- AHC_real[-small.index]
    }

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
            Q_ind <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[5]];
            Q_ind_p <- 1-pchisq(Q_ind, df = 1);

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
            AHC_new <- AHC.IV(betaX[real_snp], betaY[real_snp], seX[real_snp], seY[real_snp],
                              n, tuning, weak);

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
        for (k_iter in 1:length(AHC_real)){
          leaf <- AHC_real[[k_iter]];
          Q_ind <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[5]];
          Q_ind_p <- 1-pchisq(Q_ind, df = 1);

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
          AHC_new <- AHC.IV(betaX[real_snp], betaY[real_snp], seX[real_snp], seY[real_snp],
                            n, tuning, weak);

          AHC_real <- AHC_new ## update AHC_real;
          for (k_new in 1:length(AHC_new)){
            AHC_real[[k_new]] <- real_snp[AHC_new[[k_new]]];
          }

          junk_new <- setdiff(real_snp, sort(unlist(AHC_real)));
          AHC_junk <- sort(c(AHC_junk, junk_new));

          # trimming the small clusters
          small.index = which(lapply(AHC_real, length) <= smallcluster)
          if (length(small.index) != 0){
            AHC_junk_new <- sort(unlist(AHC_real[small.index]))
            AHC_junk = sort(union(AHC_junk_new, AHC_junk))
            AHC_real <- AHC_real[-small.index]
          }
        }

      }

    }
  }

  Nr_real <- length(AHC_real); # number of substantive clusters
  
  # clustering result, cluster id for junk SNPs = 0
  AHC_cluster <- rep(0, pz);
  
  for (k in 1:Nr_real){
    snp_index <- AHC_real[[k]];
    AHC_cluster[snp_index] <- k;
  }
  
  Nr_all <- length(unique(AHC_cluster)); # number of total clusters including the junk cluster.

  # Post-clustering estimation
  AHC_results <- matrix(NA, nrow = Nr_real, ncol = 7);
  colnames(AHC_results) <- c("ID", "length", "beta", "se", "t", "Qp", "I^2");

  for (k_r in 1:Nr_real){
    AHC_results[k_r, "ID"] <- k_r;
    AHC_results[k_r, "length"] <- length(AHC_real[[k_r]]);

    leaf_kr <- AHC_real[[k_r]];
    results_kr <- Qsq(betaX[leaf_kr], betaY[leaf_kr], seX[leaf_kr], seY[leaf_kr]);
    AHC_results[k_r, "beta"] <- results_kr[[2]];
    AHC_results[k_r, "se"] <- results_kr[[3]];
    AHC_results[k_r, "t"] <- abs(AHC_results[k_r, "beta"]/AHC_results[k_r, "se"]);
    AHC_results[k_r, "Qp"] <- pchisq(results_kr[[1]], df = length(leaf_kr) - P,
                                     lower.tail = FALSE);
    AHC_results[k_r, "I^2"] <- results_kr[[4]]
  }

  confidence_interval <- matrix(NA, nrow = Nr_real, ncol = 2);
  rownames(confidence_interval) <- paste0("Cluster", 1:Nr_real);
  colnames(confidence_interval) <- c("lower", "upper");
  confidence_interval[, "lower"] <-  AHC_results[, "beta"]-abs(qnorm(alpha/2))*AHC_results[, "se"];
  confidence_interval[, "upper"] <-  AHC_results[, "beta"]+abs(qnorm(alpha/2))*AHC_results[, "se"];

  # The null clusters
  null_index <- which(pchisq((AHC_results[,"t"])^2, df = 1, lower.tail = FALSE) > tuning) # Clusters that do not give significant estimates are classified as the null cluster.

  # Report results
  results = list( Cluster_number = Nr_all,
                  Cluster_number_real = Nr_real,
                  AHC_cluster = AHC_cluster,
                  Null_cluster = null_index,
                  Junk_cluster = AHC_junk,
                  F = F_bar,
                  AHC_results = AHC_results,
                  confidence_interval = confidence_interval)
  return(results)

}

# Internal Function (1), not to be called by users.
AHC.IV <- function(betaX, betaY, seX, seY, n, tuning = 0.1/log(n), weak = TRUE){

  # Define Constants
  pz <- length(betaY);
  P <- 1; ## this is for single outcome--single exposure case

  # Functions for the Q Test
  if (weak){
    # weak-IV robust Q
    Qsq = function(bx,by,sx,sy){
      biv        = by/bx
      k          = length(bx)
      df         = k - 1
      PL2 = function(a){
        b = a[1]
        w = 1/(sy^2/bx^2 + (b^2)*sx^2/bx^2)
        q = sum(w*(biv - b)^2)
      }
      Bhat = optimize(PL2,interval=c(-10,10))$minimum
      W = 1/(sy^2/bx^2 + (Bhat^2)*sx^2/bx^2)
      QE = sum(W*(biv - Bhat)^2)
      QE_ind = W*(biv - Bhat)^2
      QEp = 1-pchisq(QE,df)
      Isq = (QE - df)/QE
      Isq =  max(0,Isq)
      sehat = sqrt(solve(sum(W)))
      return(list(QE, Bhat, sehat, Isq, QE_ind))
    }
  }else{
    Qsq = function(bx,by,sx,sy){
      k          = length(bx)
      df         = k - 1
      y          = matrix(by/bx, nrow = k)
      s          = matrix(sy/bx, nrow = k)
      w          = 1/s^2; sum.w  = sum(w)
      mu.hat     = sum(y*w)/sum.w
      Q          = sum(w*(y-mu.hat)^2)
      se.hat     = sqrt(solve(sum.w))
      Isq        = (Q - df)/Q
      Isq        =  max(0,Isq)
      return(list(Q, mu.hat, se.hat, Isq))
    }
  }

  # AHC procedure
  comb <- combn(pz, P, repeats.allowed=F); Len <- ncol(comb);
  bvec <- matrix(betaY/betaX, nrow = pz);
  sevec <- matrix(abs(seY/betaX), nrow = pz);

  dend = hierclust(matrix(bvec, ncol = 1), matrix(1/sevec^2, ncol = 1));

  # This is to generate the selection path, i.e. the "tree" (algorithm 1)
  trees = list();
  for (b in 1:(Len - 1)){
    dstep <- cutree(dend, k=b);
    treesb = list();
    for (bb in 1:length(unique(dstep))){
      treesb[[bb]] = comb[, which(dstep == (unique(dstep))[bb])];
    }
    trees[[b]] = treesb;
  }

  # Downward testing
  start <- 2;
  clusterIV = clusternames <- list();
  while (start <= length(trees)) {
    branch = trees[[start]];
    index_single = which(lapply(branch, length) == 1);
    if (length(index_single) != 0) branch = branch[-index_single];
    for (k in 1:length(branch)){
      leaf <- sort(unique(as.vector(branch[[k]])));
      check.subset <- function(y){all(leaf %in% y)};
      index = any(unlist(lapply(clusterIV, check.subset)));
      if (length(leaf)< pz && index == FALSE){
        Q.branch <- Qsq(betaX[leaf], betaY[leaf], seX[leaf], seY[leaf])[[1]];
        p.branch <- pchisq(Q.branch, df = (length(leaf) - P), lower.tail = FALSE);
        if (p.branch > tuning){
          clusterIV[[length(clusterIV)+1]] <- leaf;
          clusternames[[length(clusternames)+1]] <- paste0('cluster', length(clusternames)+1)
        }
      }
    }
    clusterIV;
    start <- start + 1;
  }

  return(clusterIV)
}

# Internal Function (2): the dissimilarity matrix with weighting. Not to be called by users.
dissim <- function(a, wt) {
  # Inputs.   a: matrix, for which we want distances on rows,
  #           wt: masses of each row.
  # Returns.  matrix of dims. nrow(a) x nrow(a) with wtd. sqd. Eucl. distances.
  # FM, 2003/11/16

  n <- nrow(a)
  m <- ncol(a)
  adiss <- matrix(0, n, n)

  for (i1 in 2:n) {
    adiss[i1,i1] <- 0.0
    for (i2 in 1:(i1-1)) {
      adiss[i1,i2] <- 0.0
      for (j in 1:m) {
        # We use the squared Euclidean distance, weighted.
        adiss[i1,i2] <- adiss[i1,i2] + (wt[i1]*wt[i2])/(wt[i1]+wt[i2]) *
          (a[i1,j]-a[i2,j])^2
      }
      adiss[i2,i1] <- adiss[i1,i2]
    }
  }
  adiss
}

# Internal Function (3), not to be called by users.
getnns <- function(diss, flag) {
  # Inputs.  diss: full distance matrix.
  #          flag: "live" rows indicated by 1 are to be processed.
  # Returns. List of: nn, nndiss.
  #          nn:   list of nearest neighbor of each row.
  #          nndiss: nearest neigbbor distance of each row.
  # FM, 2003/11/16

  nn <- rep(0, nrow(diss))
  nndiss <- rep(0.0, nrow(diss))
  MAXVAL <- 1.0e12
  if (nrow(diss) != ncol(diss)) stop("Invalid input first parameter.")
  if (nrow(diss) != length(flag)) stop("Invalid inputs 1st/2nd parameters.")
  # if (nrow(diss) != length(nn)) stop("Invalid inputs 1st/3rd parameters.")
  # if (nrow(diss) != length(nndiss)) stop("Invalid inputs 1st/4th parameters.")

  for (i1 in 1:nrow(diss)) {
    if (flag[i1] == 1) {
      minobs <- -1
      mindis <- MAXVAL
      for (i2 in 1:ncol(diss)) {
        if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
          mindis <- diss[i1,i2]
          minobs <- i2
        }
      }
      nn[i1] <- minobs
      nndiss[i1] <- mindis
    }
  }
  list(nn = nn, nndiss = nndiss)
}

# Internal Function (3): generate the clustering path. Not to be called by users.
hierclust <- function(a, wt) {

  MAXVAL <- 1.0e12

  n <- nrow(a)
  diss <- dissim(a, wt)                      # call to function dissim
  flag <- rep(1, n)                          # active/dead indicator
  a <- rep(0, n-1)                           # left subnode on clustering
  b <- rep(0, n-1)                           # right subnode on clustering
  ia <- rep(0, n-1)                          # R-compatible version of a
  ib <- rep(0, n-1)                          # R-compatible version of b
  lev <- rep(0, n-1)                         # level or criterion values
  card <- rep(1, n)                          # cardinalities
  mass <- wt
  order <- rep(0, n)                         # R-compatible order for plotting

  nnsnnsdiss <- getnns(diss, flag)           # call to function getnns
  clusmat <- matrix(0, n, n)                 # cluster memberships
  for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition

  for (ncl in (n-1):1) {                      # main loop
    # check for agglomerable pair
    minobs <- -1;
    mindis <- MAXVAL;
    for (i in 1:n) {
      if (flag[i] == 1) {
        if (nnsnnsdiss$nndiss[i] < mindis) {
          mindis <- nnsnnsdiss$nndiss[i]
          minobs <- i
        }
      }
    }
    # find agglomerands clus1 and clus2, with former < latter
    if (minobs < nnsnnsdiss$nn[minobs]) {
      clus1 <- minobs
      clus2 <- nnsnnsdiss$nn[minobs]
    }
    if (minobs > nnsnnsdiss$nn[minobs]) {
      clus2 <- minobs
      clus1 <- nnsnnsdiss$nn[minobs]
    }
    # So, agglomeration of pair clus1 < clus2 defines cluster ncl

    #------------------------------------ Block for subnode labels
    a[ncl] <- clus1                       # aine, or left child node
    b[ncl] <- clus2                       # benjamin, or right child node
    # Now build up ia, ib as version of a, b which is R-compliant
    if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
    if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
    if (card[clus1] > 1) {                # left child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
      }
      ia[ncl] <- n - lastind             # label of non-singleton
    }
    if (card[clus2] > 1) {                # right child is non-singleton
      lastind <- 0
      for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
        if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
      }
      ib[ncl] <- n - lastind             # label of non-singleton
    }
    if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
      left <- min(ia[ncl],ib[ncl])
      right <- max(ia[ncl],ib[ncl])
      ia[ncl] <- left                    # Just get left < right
      ib[ncl] <- right
    }
    #--------------------------------------------------------------------

    lev[ncl] <- mindis
    for (i in 1:n) {
      clusmat[i,ncl] <- clusmat[i,ncl+1]
      if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
    }
    # Next we need to update diss array
    for (i in 1:n) {
      if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
        diss[clus1,i] <-
          ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,i] +
          ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus2,i] -
          (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,clus2]
        diss[i,clus1] <- diss[clus1,i]
      }
    }
    mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
    card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
    # Cluster label clus2 is knocked out; following not nec. but no harm
    flag[clus2] <- 0
    nnsnnsdiss$nndiss[clus2] <- MAXVAL
    mass[clus2] <- 0.0
    for (i in 1:n) {
      diss[clus2,i] <- MAXVAL
      diss[i,clus2] <- diss[clus2,i]
    }
    # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
    # i.e. nearest neighbors and the nearest neigh. dissimilarity
    nnsnnsdiss <- getnns(diss, flag)
  }

  temp <- cbind(a,b)
  merge2 <- temp[nrow(temp):1, ]
  temp <- cbind(ia,ib)
  merge <- temp[nrow(temp):1,]
  dimnames(merge) <- NULL
  # merge is R-compliant; later suppress merge2

  #-------------------------------- Build R-compatible order from ia, ib
  orderlist <- c(merge[n-1,1], merge[n-1,2])
  norderlist <- 2
  for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
    for (i2 in 1:norderlist) {       # Scan orderlist
      if (orderlist[i2] > 0) {     # Non-singleton to be expanded
        tobeexp <- orderlist[i2]
        if (i2 == 1) {
          orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[2:norderlist])
        }
        if (i2 == norderlist) {
          orderlist <- c(orderlist[1:(norderlist-1)],
                         merge[tobeexp,1],merge[tobeexp,2])
        }
        if (i2 > 1 && i2 < norderlist) {
          orderlist <- c(orderlist[1:(i2-1)],
                         merge[tobeexp,1],merge[tobeexp,2],
                         orderlist[(i2+1):norderlist])
        }
        norderlist <- length(orderlist)
      }
    }
  }
  orderlist <- (-orderlist)
  class(orderlist) <- "integer"

  xcall <- "hierclust(a,wt)"
  class(xcall) <- "call"
  #clusmat=clusmat
  #labels=as.character(1:n)

  retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
                  labels=dimnames(a)[[1]],method="minvar",call=xcall,
                  dist.method="euclidean-factor")
  retlist <- list(merge=merge,height=lev[(n-1):1],order=orderlist)
  class(retlist) <- "hclust"
  retlist
}
