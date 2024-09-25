
#' Detecting Changes Jointly in the Eigenvalues of the Covariance Operator of the Functional Data
#'
#' This function tests and detects changes jointly in the eigenvalue of
#'  the covariance operator.
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param d Number of eigenvalues to include in testing.
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter for the estimation of the long run covariance matrix. The default
#' value is \code{h=2}.
#' @param mean_change If \code{TRUE} then the data is centered considering the change in the mean function.
#' @param delta Trimming parameter to estimate the covariance function using partial sum estimates.
#'
#' @export
#' @return
#' \item{\code{pvalue}}{
#'  Approximate p value for testing whether there is a significant change in the desired eigenvalue of the covariance operator
#' }
#' \item{\code{change}}{
#'  Estimated change location
#' }
#' \item{\code{eval_before}}{
#'  Estimated eigenvalues before the change
#' }
#' \item{\code{eval_after}}{
#'  Estimated eigenvalues after the change
#' }
#' @details This function dates and detects changes in the joint eigenvalues
#'  that is defined by \code{d} of the covariance function. The critical
#'  values are approximated via \code{M} Monte Carlo simulations.
#'
#' @references Aue, A, G Rice, and O S\"{o}nmez. “Structural Break Analysis
#'  for Spectrum and Trace of Covariance Operators.” Environmetrics
#'  (London, Ont.) 31, no. 1 (2020). https://doi.org/10.1002/env.2617.
#'
#' @examples
#' joint_eigen_change(
#'   generate_brownian_motion(200,v=seq(0,1,length.out=20)), 1)
#' joint_eigen_change(funts(electricity), 2)
joint_eigen_change <- function(X, d, h =2, CPs = NULL,
                               delta = 0.1, M = 1000,
                               K=bartlett_kernel){
  X <- .check_data(X)
  X <- center(X, CPs=CPs)

  N <- ncol(X$data)
  D <- nrow(X$data)

  Cov_op <- .partial_cov(X, 1)

  eig2 <- array(dim = c(D,D,D))
  for(j in 1:D){
    eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
  }
  thetas <- matrix(nrow=N,ncol=D)
  X2 <- array(dim = c(D,D,N))
  for(i in 1:N){
    X2[,,i] <- X$data[,i] %*% t(X$data[,i]) - Cov_op$coef_matrix

    for(j in 1:D){
      thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
    }
  }
  Sigma_d <- .long_run_var(t(thetas[,1:d]),h=h,K = K)

  asymp <- function(N, delta, d){
    s <- floor(N * delta)
    WW <- matrix(0, d, N-s)
    for(i in 1:d){
      WW[i, ] <- sde::BBridge(0, 0, 0, 1, N)[(s+1):N]
    }

    max(colSums(WW^2))
  }
  Values_j <- sapply(1:M, function(i,N,delta,d) asymp(N, delta,d),
                     N=N, delta=delta, d=d)

  s <- floor(delta*N)
  Tn_j <-  c(rep(0,s))
  for(k in (s+1):N){
    lam_i <- .partial_cov(X, k/N)$eigen_val[1:d]
    kapa <- lam_i - (k/N) * Cov_op$eigen_val[1:d]
    Tn_j[k] <- N * t(kapa) %*% solve(Sigma_d) %*% kapa
  }
  Sn_j <- max(Tn_j)
  k_star <- min(which.max(Tn_j))
  p_j <- sum(Sn_j <= Values_j) / M # Compute p-value

  l1 <- pca(X$data[,1:k_star])$sdev[1:d]
  l2 <- pca(X$data[,(1+k_star):N])$sdev[1:d]

  list(change = k_star,
       pvalue = p_j,
       eval_before = l1,
       eval_after = l2)
}
