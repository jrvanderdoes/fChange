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
#' #specific_eigen_change(generate_brownian_motion(500),1)
#' #specific_eigen_change(funts(electricity),1)
#' #specific_eigen_change(funts(electricity),2)
specific_eigen_change <- function(X, component, h=2, CPs = NULL,
                                  delta = 0.1, M = 1000,
                                  K=bartlett_kernel){
  stop('Errors still remaining')
  set.seed(1234)
  fdata <- fun_IID(n=200, nbasis=21)
  X <- fda::eval.fd(seq(0,1,length.out=21),fdata)
  component <- 2
  h <- 2
  CPs <- NULL
  delta = 0.1
  M = 1000
  K=bartlett_kernel

  # X <- electricity
  #TODO:: Size seems very off..
  X <- .check_data(X)
  X <- center(X, CPs = CPs)

  N <- ncol(X$data)
  D <- nrow(X$data)

  # X$data is the evaled cdata
  #   Get autocov at lag 0, and its eigenval/functions
  Cov_op <- .partial_cov(X, 1)
  PSI <- Cov_op$coef_matrix
  Phi <- Cov_op$eigen_fun
  lambda <- Cov_op$eigen_val

  # Cov_op <- pca(X,TVE=1)
  # PSI <- Cov_op$x
  # Phi <- Cov_op$sdev
  # Phi <- Cov_op$rotation

  Projections <- Phi %*% X$data
  Proj_sq <- Projections^2
  Psi_diag <- diag(PSI)
  THETA <- Proj_sq - Psi_diag
  theta <- matrix(THETA[component,], ncol = N, nrow = 1)
  Sigma_d <- .long_run_var(theta, h=2, K=K)

  Values <- sapply(1:M, function(k)
    max(sde::BBridge(0, 0, 0, 1, N)[(floor(delta*N)+1):N]^2))

  s <- floor(delta*N)
  Tn <- c(rep(0,s))
  for (k in (s+1):N){
    .partial_cov1(X, k/N)$eigen_val[component]-
    lam_i <- .partial_cov1(X, k/N)$eigen_val[component]
    Tn[k] <- ( N*( lam_i - (k/N)*lambda[component] )^2 ) / Sigma_d
  }
  Sn <- max(Tn)
  # p <- ecdf(Values)(Sn)
  p <- sum(Sn <= Values) / M # Compute p-value
  p
  Sn
  Values[1:20]
  k_star <- min(which.max(Tn))

  list(change = k_star, pvalue = p)
}
