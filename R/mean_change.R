
#' Fully Functional Mean Change Point Analysis
#'
#' This function tests whether there is a significant change in the mean
#'     function of the functional data, and it will give an estimate for the
#'     location of the change. The procedure is based on the standard L-2 norm
#'     and hence does not depend on any dimension reduction technique such as
#'     fPCA.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param M (Optional) Number of Monte Carlo simulations used to get the critical
#'     values. The default value is \code{M=1000}
#' @param  h (Optional) The window parameter parameter for the estimation of the
#'     long run covariance kernel. The default value is \code{h=0}, i.e., it
#'     assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#' @param alpha (Optional) Numeric value indicating significance. Defaults to 0.05
#' @param inc.pval (Optional) Boolean to indicate if pval should also be returned
#' @param ... Unused, just for use in other functions
#'
#' @export
#' @return If inc.pval is false, Numeric of CP location or NA, otherwise
#'     location or NA and the p-value
#'
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#'
#' @examples
#' mean_change(electricity,M=250)
mean_change <- function(data, M=1000, h=0, K = bartlett_kernel, alpha=0.05,
                     inc.pval=FALSE, ...){
  n <- ncol(data)
  Sn2 <- rep(0,n)

  # CUSUM
  for(k in 2:n){
    Sn2[k]= sum((rowSums(data[,1:k]) - (k/n)*rowSums(data))^2)
  }
  Sn2 <- Sn2/n

  Tn <- max(Sn2)
  k.star <- min(which(Sn2==max(Sn2)))

  ## Estimate eigenvalues (lambda_i, 1<=i<=d)
  Ceps <- .estimateCeps(data, h, K)
  lambda <- eigen(Ceps/n)$values

  values_sim <- sapply(1:M, function(k, lambda, n) .asymp_dist(n, lambda),
                       lambda=lambda, n=n)
  p <- sum(Tn <= values_sim)/M # Compute p-value

  # Just return CP for now
  if(p<=alpha){
    return_val <- k.star
  }else{
    return_val <- NA
  }

  if(inc.pval){
    return_val <- c(return_val, p)
  }

  return_val
}


#' Asymptotic Distribution
#'
#' Define and simulate asymptotic distribution based on Brownian bridge and
#'     eigenvalues.
#'
#' @param n Integer indicating the number of observations for the Brownian bridge
#' @param lambda Vector of numerics indicating eigenvalues
#'
#' @return Numeric indicating max value of the Brownian bridge
#'
#' @noRd
#'
#' @examples
#' .asymp_dist(200,1:5)
#' .asymp_dist(200,1:5)
.asymp_dist <- function(n, lambda){
  BridgeLam <- matrix(0,length(lambda),n)
  for(j in (1:length(lambda))){
    BridgeLam[j,] <- lambda[j]*(sde::BBridge(x=0,y=0,t0=0,T=1,N=n-1)^2)
  }
  max(colSums(BridgeLam))
}


#' Estimate Long-run Covariance Kernel
#'
#' This (internal) function estimates the long-run covariance kernel. That is,
#'     \eqn{C_{\epsilon}(t,t') = \sum_{l=-\inf}^{\inf} \text{Cov}(\epsilon_0(t),
#'     \epsilon_l(t'))} with error sequence \eqn{(\epsilon_i : i \in \mathbb{Z})}.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param  h The window parameter parameter for the estimation of the long run covariance kernel. The default
#' value is \code{h=0}, i.e., it assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#'
#' @return Data.frame of numerics with dim of ncol(data) x ncol(data), that is
#'     symmetric.
#'
#' @noRd
#'
#' @examples
#' # This is an internal function, see use in mean_change.
.estimateCeps <- function(data, h, K){
  N <- ncol(data)
  D <- nrow(data)
  Ceps <- matrix(NA,nrow=D, ncol=D)

  data <- .centerData(data)

  for (k in 1:D) {
    for (r in k:D) {
      # Multiple all observations taken at the same point in time across FDs
      s <- as.numeric(data[k,]) %*% as.numeric(data[r,])
      if (h > 0) {
        for (i in 1:h) {
          # Don't fully understand
          a <- as.numeric(data[k, 1:(N - i)]) %*% as.numeric(data[r, (i + 1):N])
          a <- a + as.numeric(data[r, 1:(N - i)]) %*% as.numeric(data[k, (i + 1):N])
          s <- s + K(i, h) * a
        }
      }
      Ceps[k, r] <- Ceps[r, k] <- s
    }
  }

  Ceps
}


#' Center Data
#'
#' This (internal) function centers each time point.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#'
#' @return Numeric data.frame of the data, but with centered values at each time point
#'
#' @noRd
#'
#' @examples
#' # This is an internal function, see use in .estimateCeps.
.centerData <- function(data){
  data - rowMeans(data)
}


#' Compute CUSUM Statistic for Mean Change
#'
#' This function is used to compute the CUSUM statistic for a mean change of
#'     FD observations.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param k Numeric indicating the change point of interest
#' @param ... Unused, just for use in other functions
#'
#' @return Numeric indicating the CUSUM value at the given K value
#' @export
#'
#' @examples
#' \dontrun{
#' # Null Example
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0,0),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#'
#' compute_mean_stat(data_KL, 100)
#'
#' # Mean CP Example
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0,0.2),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#'
#' compute_mean_stat(data_KL, 100)
#' }
compute_mean_stat <- function(data, k, ...){
  n <- ncol(data)

  sum((rowSums(as.data.frame(data[,1:k])) - (k/n)*rowSums(as.data.frame(data)))^2)/n
}


#' Compute CUSUM Mean Statistic Cutoff
#'
#' This function computes the cutoff for CUSUM statistic of a mean change
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param alpha Numeric value indicating significance
#' @param  h (Optional) The window parameter parameter for the estimation of the
#'     long run covariance kernel. The default value is \code{h=0}, i.e., it
#'     assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#' @param M (Optional) Number of Monte Carlo simulations used to get the critical
#'     values. The default value is \code{M=1000}
#' @param ... Unused, just for use in other functions
#'
#' @return Numeric indicating the cutoff value for CUSUM mean statistic
#' @export
#'
#' @examples
#' \dontrun{
#' # Null Example
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0,0),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#'
#' compute_mean_cutoff(data_KL, 0.05)
#'
#' # Mean CP Example
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0,0.2),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#'
#' compute_mean_cutoff(data_KL, 0.05)
#' }
compute_mean_cutoff <- function(data, alpha, h=0, K=bartlett_kernel,
                                M=1000, ...){
  n <- ncol(data)

  ## Estimate eigenvalues (lambda_i, 1<=i<=d)
  Ceps <- .estimateCeps(data, h, K)
  lambda <- eigen(Ceps/n)$values

  values_sim <- sapply(1:M, function(k, lambda, n) .asymp_dist(n, lambda),
                       lambda=lambda,n=n)

  as.numeric(stats::quantile(values_sim,1-alpha))
}
