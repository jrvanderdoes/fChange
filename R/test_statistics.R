
#' Compute Tn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Tn = \int_0^1 \int |Z_n(v,x)|^2 dQ(v) dx}
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD.
#' @param k (Optional) Integer indicating if Defaults to NULL, but can be an
#' integer in \eqn{[1, ncol(X)]}.
#'
#' Indicates the candidate change point to investigate (uses `compute_Mn()`
#'  for selection). If k=NULL then get the test statistic for entire sample.
#'  Default is NULL.
#'
#' @param M (Optional) Numeric indicating the number of vectors used to span W.
#'  Defaults to 10,000.
#' @param W (Optional) Data.frame of numerics with rows for evaluated values
#'     and columns indicating function that contains the vectors to span the
#'     space. Note the number of columns should match that of X. Default is NULL,
#'     which generates the vectors in function. This overrides space and M if
#'     not NULL.
#' @param space (Optional) String indicating the space to integrate against.
#'  Default is "BM", or Brownian motion. See `computeSpaceMeasuringVectors()`
#'  for more information.
#' @param ... Parameters are not passed anywhere. Just added for use in other
#'  functions.
#'
#' @return Numeric value for the test statistic of entire sample or candidate
#'     change point.
#' @export
#'
#' @examples
#' compute_Tn(electricity)
compute_Tn <- function(X, M=100000, W=NULL, space='BM', ...){
  n <- ncol(X)

  if(is.null(W)){
    W <- computeSpaceMeasuringVectors(M = M, X = X, space=space)
  }else{
    M <- ncol(W)
  }

  intVal <- sapply(1:M, function(v,W,X1,n){
    .approx_int( (abs(.Zn(W[,v],X1)))^2 ) #RH Int
  },W=W,X1=X,n=n)

  1/M * sum(intVal)
}


#' Compute Mn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Mn = sup_{x\in (0, 1)} \int |Z_n(v,x)|^2 dQ(v)}
#'
#' @inheritParams compute_Tn
#' @param which.Mn (Optional) Boolean which indicates if the location or test
#'  statistic is of interest.
#'
#' @return A numeric value for the test statistic for entire sample or candidate
#'     change point.
#' @export
#'
#' @examples
#' compute_Mn(electricity)
compute_Mn <- function(X, M=10000, W=NULL, space='BM', ...){
  n <- ncol(X)

  if(is.null(W)){
    W <- computeSpaceMeasuringVectors(M = M, X = X, space=space)
  }else{
    M <- ncol(W)
  }

  intVal <- sapply(1:M, function(v,W,X1,n){
    (abs(.Zn(W[,v],X1)))^2
  },W=W,X1=X,n=n)
  return_value <- 1/M * rowSums(intVal)

  list('value'=max(return_value),
       'location'=which.max(return_value),
       'allValues'=return_value)
}


#' Compute Zn Statistic
#'
#' This computes Zn. See usage in compute_Tn.
#'
#' @param v Vector of numerics with length of each curve in X, nrow(X).
#' @param X Data.frame of functional data. Observations are columns, resolution
#'  is the rows.
#'
#' @return Numeric value for Zn, where \eqn{Zn(v,x)=\sqrt{n} (\hat{f}_n(v,x) -
#' \frac{\lfloor n x \rfloor}{n} \hat{f}_n(v,1))}.
#'
#' @noRd
.Zn <- function(v,X){
  n <- ncol(X)
  fhat_vals <- .fhat_all(X,v)
  #fhat_full <- sum(fhat_vals)

  sqrt(n) * ( fhat_vals -
               1:n/n * fhat_vals[length(fhat_vals)] )
}


#' Get fhat Estimates
#'
#' This computes the values in fHat, i.e. cutoff for each possible value.
#'  Note, a starting 0 is omitted. See usage in .Zn.
#'
#' @inheritParams .Zn
#'
#' @return Vector of numerics for \eqn{\hat{f}(v,x)}
#'
#' @noRd
.fhat_all <- function(X,v){
  cumsum(exp(1/nrow(X)*complex(imaginary = 1)* t(X) %*% v)) / ncol(X)
}
