

#' Testing changes in the trace of the covariance operator in functional data
#'
#' This function tests and detects changes in the trace of the covariance operator.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param CPs If not NULL then the data is centered considering the changes.
#' @param delta Trimming parameter to estimate the covariance function using partial sum estimates.
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#'  value is \code{h=2}.
#'
#'
#'
#'
#' @return
#' \item{\code{pvalue}}{
#'  Approximate p value for testing whether there is a significant change in the desired eigenvalue of the covariance operator
#' }
#' \item{\code{change}}{
#'  Estimated change location
#' }
#' \item{\code{trace_before}}{
#'  Estimated trace before the change
#' }
#' \item{\code{trace_after}}{
#'  Estimated trace after the change
#' }
#' @export
#'
#' @details This function dates and detects changes in trace of the covariance
#'  function. This can be interpreted as the changes in the total variation
#'  of the the functional data. Trace is defined as the infinite sum of
#'  the eigenvalues of the covariance operator and for the sake of
#'  implementation purpose, the sum is truncated up to the total number
#'  of basis functions that defines the functional data at hand. The critical
#'  values are approximated via \code{M} Monte Carlo simulations.
#'
#' @references Aue, A, G Rice, and O Sönmez. “Structural Break Analysis
#'  for Spectrum and Trace of Covariance Operators.” Environmetrics
#'  (London, Ont.) 31, no. 1 (2020). \url{https://doi.org/10.1002/env.2617.}
#'
#' @examples
#' trace_change(generate_brownian_motion(200,v=seq(0,1,0.05)))
#' trace_change(electricity)
trace_change <- function (X, CPs = NULL, M = 1000){
  X <- .check_data(X)
  X <- center(X, CPs=CPs)

  N <- ncol(X$data)
  D <- nrow(X$data)

  Cov_op <- .partial_cov(X, 1)
  lambda <- Cov_op$eigen_val
  T_1 <- sum(lambda)
  Xi <- sapply(1:N, function(i) X$data[,i] %*% X$data[,i])
  sigma_sq <- sandwich::lrvar(Xi, prewhite = F)
  sigma <- sqrt(sigma_sq)
  Values <- sapply(1:M, function(k) max(abs(sde::BBridge(0, 0, 0, 1, N-1))))

  Tn <- rep(0,N)
  for (k in 1:N) {
    T_x <- sum(Xi[1:k])/N
    Tn[k] <- (1/sigma) * abs(T_x - (k)/N * T_1)
  }
  Sn <- max(Tn)
  # p <- ecdf(Values)(Sn)
  p <- sum(Sn <= Values) / M # Compute p-value
  k_star <- min(which.max(Tn))

  data1_pca <- pca(funts(X$data[,1:k_star]),TVE=1)
  data2_pca <- pca(funts(X$data[,(k_star+1):N]),TVE=1)
  tr_before <- sum(data1_pca$sdev)
  tr_after <- sum(data2_pca$sdev)
  list(change = k_star/N,
       location = k_star,
       pvalue = p,
       trace_before = tr_before,
       trace_after = tr_after)
}


#' Partial Sample Estimates of the Covariance Function in Functional Data Analysis
#'
#' This function computes the partial sum estimate of the covariance function
#'  in functional data analysis. It also generate the eigenvalues and
#'  eigenfunctions of the parial sum estimate, along with the coefficient
#'  matrix of the estimated covariance function.
#'
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param x Fraction of the sample size that the partial sum is computed.
#'  This input must be in (0,1]. The default is \code{x=1} which corresponds
#'  to the regular covariance function estimate in functional data analysis.
#'
#' @export
#'
#' @return
#'\item{\code{eigen_val}}{
#' Eigenvalues of the partial sum estimate of the covariance function
#'}
#'\item{\code{eigen_fun}}{
#' Eigenfunctions of the partial sum estimate of the covariance function
#'}
#'\item{\code{coef_matrix}}{
#' Coefficient matrix of the partial sum estimate of the covariance function
#'}
#'
#'@details This function simply estimates the covariance function based on the partial sum of the centered functional observations.
#'The length of the sum is determined by \code{x}, and when \code{x=1} this estimate corresponds to the regular covariance function
#'estimates using the whole sample.
#'
#' @seealso \code{\link{pca.fd}}
#'
#' @examples
#' .partial_cov(electricity)
#' # Estimated eigenvalues
#' e1 = eigen(autocov_approx_h(electricity,0))
#' e2 = .partial_cov(funts(electricity),x=1)
#' sum(round(e1$values,4) != round(e2$eigen_val,4) )
#' sum(round(e1$vectors,4) != round(e2$eigen_fun,4) )
#' sum(round(autocov_approx_h(electricity,0),4) != round(e2$coef_matrix,4) )
#' # e1 and e2 will both estimate the eigenvalues of the covariance
#' # operator based on the whole sample
#'
#' # estimates using only 90% of the data
#' Cov = .partial_cov(electricity, 0.9)
.partial_cov <- function(X, x = NULL){
  X <- .check_data(X)
  if (is.null(x)) x <- 1
  if (x>1 | x<=0) stop("x should be in (0,1]")

  n <- ncol(X$data)
  D <- nrow(X$data)
  k <- floor(n*x)

  X <- center(X)
  C <- matrix(0, D, D)
  for (i in 1:k){
    C <- C + X$data[,i] %*% t(X$data[,i])
  }
  E <- eigen(C/n)
  e_val <- E$values
  e_fun <- E$vectors # fd(E$vectors, fdobj$basis)

  list(eigen_val = E$values,
       eigen_fun = E$vectors,
       coef_matrix = C/n)
}

