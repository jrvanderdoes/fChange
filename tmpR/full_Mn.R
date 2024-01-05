# Full Match - Mn
nTrials <- 2500
pvals <- pvals1 <- pvals2 <- rep(NA,nTrials)
set.seed(1234)
for(i in 1:nTrials){
  cat(paste0(i,', '))
  X <- .generateFunctionalIID(21,100)

  tmp <- detect_changepoint_final_Mn(X,M=5, nSims = 1000,h=0)
  pvals[i] <- tmp$pval
  pvals1[i] <- tmp$pval1
  pvals2[i] <- tmp$pval2
}

library(ggplot2)
ggplot() +
  geom_line(aes(x=1:nTrials,y=dplyr::cummean(pvals2<=0.05),
                color='None'),linetype='solid') +
  geom_line(aes(x=1:nTrials,y=dplyr::cummean(pvals<=0.05),
                color='After'),linetype='dashed') +
  geom_line(aes(x=1:nTrials,y=dplyr::cummean(pvals1<=0.05),
                color='During'),linetype='dotted') +
  geom_point(aes(x=nTrials,y=mean(pvals<=0.05),color='After')) +
  geom_point(aes(x=nTrials,y=mean(pvals1<=0.05),color='During')) +
  geom_point(aes(x=nTrials,y=mean(pvals2<=0.05),color='None')) +
  theme_bw() +
  geom_hline(aes(yintercept=0.05))

detect_changepoint_final_Mn <- function(X,
                                        M = 40, nSims = 2000,
                                        h = 3,
                                        K = bartlett_kernel,
                                        space = "BM",
                                        silent = FALSE) {
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_Mn <- compute_Mn_final(X, W)
  val_Mn1 <- compute_Mn_final1(X, W)
  val_Mn2 <- compute_Mn_final2(X, W)

  MJ <- M * 100

  sqrtMat <- .compute_sqrtMat_final(X,W,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)
    # 0 0 0 0 (M times) ... .... 1 1 1 1 (M times)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    apply(results, MARGIN = 1, max)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=100)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Mn$value),
    "pval1" = 1 - stats::ecdf(gamProcess)(val_Mn1$value),
    "pval2" = 1 - stats::ecdf(gamProcess)(val_Mn2$value),
    "gamProcess" = gamProcess,
    "value" = val_Mn$value,
    "value1" = val_Mn1$value,
    "value2" = val_Mn2$value
  )
}



compute_Mn_final <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final(W,X)
  ns <- .select_n(1:n,100)
  Zn <- Zn[ns,]
  # Integrate out W
  return_value <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  list(
    "value" = max(return_value),
    "location" = which.max(return_value),
    "allValues" = return_value
  )
}


compute_Mn_final1 <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final1(W, X)
  # Integrate out W
  return_value <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  list(
    "value" = max(return_value),
    "location" = which.max(return_value),
    "allValues" = return_value
  )
}

compute_Mn_final2 <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final(W,X)
  # Integrate out W
  return_value <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  list(
    "value" = max(return_value),
    "location" = which.max(return_value),
    "allValues" = return_value
  )
}



.generateFunctionalIID <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = TRUE)

  t(fiid)
}



detect_changepoint_final <- function(X,
                                     M = 40, nSims = 2000,
                                     h = 3,
                                     K = bartlett_kernel,
                                     space = "BM",
                                     silent = FALSE) {
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_Tn <- compute_Tn_final(X, W)
  val_Tn1 <- compute_Tn_final1(X, W)
  val_Tn2 <- compute_Tn_final2(X, W)

  MJ <- M * 100

  sqrtMat <- .compute_sqrtMat_final(X,W,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    # todo / ncol(mvnorms1)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    #apply(abs(gamVals)^2, MARGIN = 1, dot_integrate)
    apply(results, MARGIN = 1, dot_integrate)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=100)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Tn),
    "pval1" = 1 - stats::ecdf(gamProcess)(val_Tn1),
    "pval2" = 1 - stats::ecdf(gamProcess)(val_Tn2),
    "gamProcess" = gamProcess,
    "value" = val_Tn,
    "value1" = val_Tn1,
    "value2" = val_Tn2
  )
}



#############################################


# compute_Mn_final <- function(X, W) {
#   n <- ncol(X)
#   M <- ncol(W)
#
#   Zn <- abs(.Zn_final(W,X))^2
#   return_value <- unname(unlist( 1 / M * rowSums(Zn) ))
#
#   list(
#     "value" = max(return_value),
#     "location" = which.max(return_value),
#     "allValues" = return_value
#   )
# }

compute_Tn_final <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final(W, X)
  ns <- .select_n(1:n,100)
  Zn <- Zn[ns,]
  # Integrate out W
  intVal <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  # Integrate observations
  dot_integrate(intVal)
}

.Zn_final <- function(W, X) {
  n <- ncol(X)
  fhat_vals <- as.matrix(.fhat_all(X, W))

  sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])
}


compute_Tn_final1 <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final(W, X)
  # Integrate out W
  intVal <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  # Integrate observations
  dot_integrate(intVal)
}

.Zn_final1 <- function(W, X) {
  fhat_vals <- as.matrix(.fhat_all(X, W))

  n <- ncol(X)
  res <- 100
  ns <- .select_n(1:ncol(X),res)

  sqrt(n) * (fhat_vals[ns,] - (ns/n) %o% fhat_vals[nrow(fhat_vals),])
}

compute_Tn_final2 <- function(X, W) {
  n <- ncol(X)
  M <- ncol(W)

  Zn <- .Zn_final(W, X)
  # Integrate out W
  intVal <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  # Integrate observations
  dot_integrate(intVal)
}

.select_n <- function(vals,n){
  vals[round(seq(1,length(vals),length.out=n))]
}

#############################################

.compute_sqrtMat_final <- function(X,W,h,K){
  x <- seq(0,1,length.out=100)

  # Compute Gamma Matrix
  covMat <- .estimCovMat_final(X, W, x, h, K)
  tryCatch({
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  }, error = function(e){
    covMat <- round(covMat, 15)
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  })

  Rfast::mat.mult(
    Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
  )
}

#' Estimate the Covariance Matrix
#'
#' This (internal) function computes the covariance matrix of some FD data.
#'
#' @inheritParams detect_changepoint
#' @param M (Optional) Integer indicating the number of vectors used to create
#'     each value in the covariance matrix. Default is 25.
#' @param W (Optional) Data.frame of vectors for measuring space. If NULL then
#'     then a Guassian measure is generated. Default is NULL.
#'
#' @return Data.frame for covariance based on Gaussian measure and given data X
#'
#' @noRd
.estimCovMat_final <- function(X, W, x=seq(0,1,length.out=100),
                               h = 3, K = bartlett_kernel) {
  M <- ncol(W)

  D11 <- D12_tD21 <- D22 <- data.frame(matrix(ncol = M, nrow = M))
  funs <- c(cos, sin)
  for (i in 1:M) {
    for (j in i:M) {
      D11[i, j] <- D11[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[1][[1]], vp = W[, j]
        )
      D12_tD21[i, j] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
        )
      if (j > i) {
        D12_tD21[j, i] <- .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, j],
          lfunp = funs[2][[1]], vp = W[, i]
        )
      }

      D22[i, j] <- D22[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[2][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
        )
    }
  }

  tmp1 <- x %*% t(x)
  minDF <- matrix(1, nrow = length(x), ncol = length(x))
  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      minDF[i, j] <- (min(x[i], x[j]) - tmp1[i, j])
    }
  }

  fastmatrix::kronecker.prod(
    as.matrix(minDF),
    rbind(
      cbind(as.matrix(D11), as.matrix(D12_tD21)),
      cbind(as.matrix(t(D12_tD21)), as.matrix(D22))
    )
  )
}
