
#' Compute Tn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Tn = \int_0^1 \int |Z_n(v,x)|^2 dQ(v) dx}
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param k (Optional) Defaults to NULL, but can be an integer in [1,ncol(X)]
#'
#' Indicates the candidate change point to investigate. If k=NULL then get the
#'     test statistic for entire sample
#'
#' @param M (Optional) Defaults to 100. The number of functions in W
#' @param W (Optional) Optional numeric matrix with rows for evaluated values
#'     and columns indicating function
#'
#' This is the respective functions used to investigate the test statistic. The
#'     NULL value means Brownian motion is used.
#'
#' @return A numeric value for the test statistic for entire sample or candidate
#'     change point.
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(200),
#'     eigsList = list(c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0),
#'     distsArray = c('Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0))
#' compute_Tn(data_KL)
#'
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-0.5,0.5),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' compute_Tn(data_KL)
compute_Tn <- function(X, k=NULL, M=100, W=NULL, ...){
  n <- length(X[1,])
  Xn <- length(X[,1])

  if(is.null(W)){
    W <- as.data.frame(sapply(rep(0,M),sde::BM, N=Xn-1))
  }else{
    W <- as.data.frame(W)
  }

  if(!is.null(k)){
    #return_value <- 1/M * sum(sapply(W, .combZnInt, X1=X, n=n, nx=0:k))
    return_value <- 1/M * rowSums(as.data.frame(sapply(W,.combZn,nx=0:n,X1=X,n=n)))[k]
  }else{
    return_value <- 1/M * sum(as.data.frame(sapply(W, .combZnInt, X1=X, n=n, nx=0:n)))
  }

  return_value
}


#' Compute Integrated Zn Squared ($\int |Zn|^2$)
#'
#' This (internal) function computes $\int |Zn|^2$.
#'
#' @param v Vector of numeric for vector in space, length(v) = nrow(X1)
#' @param X1 Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param n Integer indicating the number of FD observations
#' @param nx Vector of integer indicating indices to consider
#'
#' @return Vector of numerics for  $\int |Zn|^2$.
#'
#' @examples
#' # This is an internal function, see use in compute_Tn for use.
.combZnInt <- function(v, X1, n, nx){
  # Rounding occasionally messes up floor(nx). Try floor(116/400*400) to see

  nfHat_vals <- exp(complex(real=0,imaginary = 1) * (t(X1) %*% v))

  # Get nFhat - first value will be 0
  nfx <- c(0,sapply(nx[-1], function(x, nfHat_vals){sum(nfHat_vals[1:x])},
                    nfHat_vals=nfHat_vals))

  # Compute abs(Zn^2), Zn= sqrt(n)*(f(v,x)-floor(nx)/n*f(v,1))
  #   sqrt(n)/n = 1/sqrt(n)
  fVals <- abs(1/sqrt(n)*(nfx-nx/n*nfx[length(nfx)]))^2

  # LH Integration
  sum((fVals[-length(fVals)])) * 1/n
}


#' Compute Mn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Mn = sup_{x\in (0,1)} \int |Z_n(v,x)|^2 dQ(v)}
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param k (Optional) Defaults to NULL, but can be an integer in [1,ncol(X)]
#'
#' Indicates the candidate change point to investigate. If k=NULL then get the
#'     test statistic for entire sample
#'
#' @param M (Optional) Defaults to 100. The number of functions in W
#' @param W (Optional) Optional numeric matrix with rows for evaluated values
#'     and columns indicating function
#'
#' This is the respective functions used to investigate the test statistic. The
#'     NULL value means Brownian motion is used.
#'
#' @return A numeric value for the test statistic for entire sample or candidate
#'     change point.
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(200),
#'     eigsList = list(c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0),
#'     distsArray = c('Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0))
#' compute_Mn(data_KL)
#'
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-0.5,0.5),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' compute_Mn(data_KL)
compute_Mn <- function(X, k=NULL, M=100, W=NULL, ...){

  n <- length(X[1,])
  xSeq <- seq(1/n, 1, 1/n)
  xN <- length(X[,1])
  if(is.null(W)){
    W <- as.data.frame(sapply(rep(0,M),sde::BM, N=xN-1))
  }


  if(!is.null(k)){
    return_value <- 1/M * rowSums(as.data.frame(sapply(W,.combZn,nx=0:n,X1=X,n=n)))[k]
  }else{
    return_value <- max(1/M * rowSums(as.data.frame(sapply(W,.combZn,nx=0:n,X1=X,n=n))))
  }

  return_value
}


#' Compute Zn Squared, |Zn|^2
#'
#' This (internal) function computes |Z_n|^2
#'
#' @param v Vector of numeric for vector in space, length(v) = nrow(X1)
#' @param X1 Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param n Integer indicating the number of FD observations
#' @param nx Vector of integer indicating indices to consider
#'
#' @return Vector of numerics for  $|Zn|^2$.
#'
#' @examples
#' # This is an internal function, see use in compute_Tn for use.
.combZn <- function(v,X1,n,nx){
  # Rounding occasionally messes up floor(nx). Try floor(116/400*400) to see

  nfHat_vals <- exp(complex(real=0,imaginary = 1) * (t(X1) %*% v))

  # Get nFhat - first value will be 0
  nfx <- c(0,sapply(nx[-1], function(x, nfHat_vals){sum(nfHat_vals[1:x])},
                    nfHat_vals=nfHat_vals))

  # Compute abs(Zn^2), Zn= sqrt(n)*(f(v,x)-floor(nx)/n*f(v,1))
  #   sqrt(n)/n = 1/sqrt(n)

  (abs(1/sqrt(n)*(nfx-nx/n*nfx[length(nfx)]))^2)[-1]
}

