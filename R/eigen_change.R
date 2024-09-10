
#' Eigenvalue Changes - Detect eigenvalue changes in the Covariance Operator
#'
#' This function tests and detects changes individual or jointly in the
#'  eigenvalue(s) of the covariance operator.
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
#' @references Aue, A, G Rice, and O Sönmez. “Structural Break Analysis
#'  for Spectrum and Trace of Covariance Operators.” Environmetrics
#'  (London, Ont.) 31, no. 1 (2020). https://doi.org/10.1002/env.2617.
#'
#' @examples
#' bm <- generate_brownian_motion(200,v=seq(0,1,length.out=20))
#' eigen_change(bm, 3, test='joint')
#' eigen_change(bm, 3, test='individual')
#' eigen_change(funts(electricity), 2, test='joint')
#' eigen_change(funts(electricity), 2, test='individual')
eigen_change <- function(X, d, h =2, CPs = NULL,
                               test = c('joint', 'individual'),
                               M = 1000,
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

  poss_tests <- c('joint', 'individual')
  test <- poss_tests[min(pmatch(test,poss_tests))]
  if(test=='joint'){
    Sigma_d <- .long_run_var(t(thetas[,1:d]),h=h,K = K)

    Values <- sapply(1:M,
                     function(i,N,d) .asymp_joint(N,d),
                     N=N, d=d)

    Tns <-  rep(0,N)
    for(k in 1:N){
      lam <- .partial_cov(X, k/N)$eigen_val[1:d]
      kapa <- lam - (k/N) * Cov_op$eigen_val[1:d]
      Tns[k] <- N * t(kapa) %*% solve(Sigma_d) %*% kapa
    }

    Sn <- max(Tns)
    k_star <- min(which.max(Tns))
    p <- sum(Sn <= Values) / M # Compute p-value

    before <- pca(X$data[,1:k_star])$sdev[1:d]
    after <- pca(X$data[,(1+k_star):N])$sdev[1:d]
  } else if(test=='individual'){
    Sigma_d <- .long_run_var(t(thetas[,d]),h=h,K = K)

    Values <- sapply(1:M, function(k)
      max(sde::BBridge(0, 0, 0, 1, N)^2))

    Tns <- rep(0,N)
    for (k in 1:N){
      lam <- .partial_cov(X, k/N)$eigen_val[d]
      Tns[k] <- ( N*( lam - (k/N)*Cov_op$eigen_val[d] )^2 ) / Sigma_d
    }

    Sn <- max(Tns)
    p <- sum(Sn <= Values) / M # Compute p-value
    k_star <- min(which.max(Tns))


    before <- pca(X$data[,1:k_star])$sdev[d]
    after <- pca(X$data[,(1+k_star):N])$sdev[d]
  }else{
    stop('test must be "joint" or "individual".', call. = FALSE)
  }

  list(change = k_star,
       pvalue = p,
       eval_before = before,
       eval_after = after)
}


#' Title
#'
#' @param N
#' @param delta
#' @param d
#'
#' @return
.asymp_joint <- function(N, d){
  WW <- matrix(0, d, N)
  for(i in 1:d){
    WW[i, ] <- sde::BBridge(0, 0, 0, 1, N-1)
  }

  max(colSums(WW^2))
}
