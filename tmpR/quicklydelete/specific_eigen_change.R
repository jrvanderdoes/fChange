#' Detecting Changes in the Eigenvalues of the Covariance Operator of the Functional Data
#'
#' This function tests and detects changes in the specific eigenvalue of
#'  the covariance operator.
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param component The eigenvalue that the componentwise test is applied to.
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter for the estimation of the long run covariance matrix. The default
#'  value is \code{h=2}.
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
#'
#' @details This function dates and detects changes in the defined eigenvalue of the covariance function. The critical values are
#'  approximated via \code{M} Monte Carlo simulations.
#'
#' @references Aue, A, G Rice, and O Sönmez. “Structural Break Analysis
#'  for Spectrum and Trace of Covariance Operators.” Environmetrics
#'  (London, Ont.) 31, no. 1 (2020). https://doi.org/10.1002/env.2617.
#'
#' @examples
#' specific_eigen_change(generate_brownian_motion(100, v=seq(0,1,length.out=20)),1)
#' specific_eigen_change(funts(electricity),1)
#' specific_eigen_change(funts(electricity),2)
specific_eigen_change <- function(X, component, h=2, CPs = NULL,
                                  delta = 0.1, M = 1000,
                                  K=bartlett_kernel){
  X <- .check_data(X)
  X <- center(X, CPs = CPs)

  N <- ncol(X$data)
  D <- nrow(X$data)

  # X$data is the evaled cdata
  #   Get autocov at lag 0, and its eigenval/functions
  Cov_op <- .partial_cov(X, 1)

  eig2 <- array(dim = c(D,D,D))
  for(j in 1:D){
    eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
  }
  thetas <- matrix(nrow=N,ncol=1)
  for(i in 1:N){
    X2 <- X$data[,i] %*% t(X$data[,i]) - Cov_op$coef_matrix

    thetas[i,1] <- sum(diag( t(eig2[,,component]) %*% X2 ))
  }
  Sigma_d <- .long_run_var(t(thetas),h=h,K = K)

  Values <- sapply(1:M, function(k)
    max(sde::BBridge(0, 0, 0, 1, N)[(floor(delta*N)+1):N]^2))

  s <- floor(delta*N)
  Tn <- c(rep(0,s))
  for (k in (s+1):N){
    lam_i <- .partial_cov(X, k/N)$eigen_val[component]
    Tn[k] <- ( N*( lam_i - (k/N)*Cov_op$eigen_val[component] )^2 ) / Sigma_d
  }
  Sn <- max(Tn)
  # p <- ecdf(Values)(Sn)
  p <- sum(Sn <= Values) / M # Compute p-value
  k_star <- min(which.max(Tn))

  list(change = k_star, pvalue = p)
}
