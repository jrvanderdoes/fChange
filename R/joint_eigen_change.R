
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
#'
#'
#' @examples
#' #joint_eigen_change(generate_brownian_motion(500), 1)
#' #joint_eigen_change(funts(electricity), 2)
joint_eigen_change <- function(X, d, h =2, CPs = NULL,
                               delta = 0.1, M = 1000,
                               K=bartlett_kernel){
  stop('Errors still remaining')
  #TODO:: Size seems very off..
  #TODO:: Refernece
  X <- .check_data(X)
  X <- center(X, CPs=CPs)

  N <- ncol(X$data)
  D <- nrow(X$data)

  Cov_op <- .partial_cov(X, 1)
  PSI <- Cov_op$coef_matrix
  Phi <- Cov_op$eigen_fun
  lambda <- Cov_op$eigen_val

  Projections <- Phi %*% X$data
  Proj_sq <- Projections^2
  Psi_diag <- diag(PSI)
  THETA <- Proj_sq - Psi_diag
  theta <- matrix(THETA[1:d,], ncol = N, nrow = d)
  Sigma_d <- .estimateCeps(theta, h=h, K=K)

  asymp <- function(){
    s <- floor(N * delta)
    WW <- matrix(0, d, N-s)
    for(i in 1:d){
      WW[i, ] <- sde::BBridge(0, 0, 0, 1, N)[(s+1):N]
    }

    max(colSums(WW^2))
  }
  Values_j <- sapply(1:M, function(i) asymp())

  s <- floor(delta*N)
  Tn_j <-  c(rep(0,s))
  for(k in (s+1):N){
    lam_i <- .partial_cov(X, k/N)$eigen_val[1:d]
    kapa <- lam_i - (k/N) * lambda[1:d]
    Tn_j[k] <- N * t(kapa) %*% solve(Sigma_d) %*% kapa
  }
  Sn_j <- max(Tn_j)
  k_star <- min(which.max(Tn_j))
  # z_j <- Sn_j <= Values_j
  # p <- ecdf(Values)(Sn)
  p_j <- sum(Sn_j <= Values_j) / M # Compute p-value
  # p_j <- length(z_j[z_j==TRUE])/length(z_j)

  l1 <- pca(X$data[,1:k_star])$sdev[1:d]
  l2 <- pca(X$data[,(1+k_star):N])$sdev[1:d]
  # l1 <- pca.fd(fdobj[1:k_star], nharm = D, centerfns = T)$values
  # l2 <- pca.fd(fdobj[(1+k_star):N], nharm = D, centerfns = T)$values

  list(change = k_star/N,
       location = k_star,
       pvalue = p_j,
       eval_before = l1,
       eval_after = l2)
}
