
#' Fully Functional Mean Change Point Analysis
#'
#' This function tests whether there is a significant change in the mean
#'     function of the functional data, and it will give an estimate for the
#'     location of the change. The procedure is based on the standard L-2 norm
#'     and hence does not depend on any dimension reduction technique such as
#'     fPCA.
#'
#' @param data Numeric data.frame with evaled points on rows and fd objects in columns
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter parameter for the estimation of the long run covariance kernel. The default
#' value is \code{h=0}, i.e., it assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#' @param alpha (Optional) Numeric value indicating significance. Defaults to 0.05
#' @param inc.pval (Optional) Boolean to indicate if pval should also be returned
#'
#' @export
#' @return If inc.pval is false, Numeric of CP location or NA, otherwise
#'     location or NA and the p-value
#'
#' @references Aue A., Rice G., Sonmez O. (2017+), \emph{Detecting and dating structural breaks in
#' functional data without dimension reduction} (https://arxiv.org/pdf/1511.04020.pdf)
#'
#' @examples
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
#'  mean_change(data_KL)
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
#' mean_change(data_KL)
mean_change <- function(data, M=1000, h=0,
                     K = bartlett_kernel, alpha=0.05,
                     inc.pval=F){
  ## Code
  n <- ncol(data)

  Sn2 <- rep(0,n)
  for(k in 2:n){
    Sn2[k]= sum((rowSums(data[,1:k]) - (k/n)*rowSums(data))^2)
  }
  Sn2 <- Sn2/n

  k.star <- min(which(Sn2==max(Sn2)))
  Tn <- max(Sn2)

  ## Use for finding change function
  #pre_CP = data[1:k.star]
  #post_CP = data[(k.star+1):N]
  #pre_mean = rowMeans(pre_CP)
  #post_mean = rowMeans(post_CP)
  #est_change_func = pre_mean - post_mean

  ## Estimate eigenvalues (lambda_i, 1<=i<=d)

  Ceps <- .estimateCeps(data, h, K)

  lambda <- eigen(Ceps/ncol(data))$values

  ## Define and simulate asymptotic distribution for p-value
  asymp_dist <- function(n, lambda){
    BridgeLam= matrix(0,length(lambda),n)
    for(j in (1:length(lambda))){
      BridgeLam[j,]=lambda[j]*(sde::BBridge(x=0,y=0,t0=0,T=1,N=n-1)^2)
    }
    max(colSums(BridgeLam))
  }

  values_sim <- sapply(1:M, function(k, lambda) asymp_dist(n, lambda),
                       lambda=lambda)
  p <- sum(Tn <= values_sim)/M

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

.estimateCeps <- function(data, h, K){
  ## Functions
  .centerData <- function(data){
    data - rowMeans(data)
  }

  ## Code
  N <- ncol(data)
  D <- nrow(data)
  data <- .centerData(data)

  Ceps <- matrix(NA,nrow=D, ncol=D)

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
