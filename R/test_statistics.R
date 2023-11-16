#' Compute Tn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Tn = \int_0^1 \int |Z_n(v,x)|^2 dQ(v) dx}
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD.
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
#' compute_Tn(electricity, M = 1000)
compute_Tn <- function(X, M = 100000, W = NULL, space = "BM", ...) {
  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }


  Zn <- .Zn(W,X, n)
  intVal <- dot_integrate_col(abs(Zn)^2)

  sum(intVal)/M
}

#' Compute Mn Test Statistic
#'
#' Function to calculate the change point test statistic for functional data
#'     \deqn{Mn = sup_{x\in (0, 1)} \int |Z_n(v,x)|^2 dQ(v)}
#'
#' @inheritParams compute_Tn
#'
#' @return A list with three elements:
#'  \itemize{
#'    \item **value**: Numeric for maximum Mn value in data.
#'    \item **location**: Numeric for location of maximum in data.
#'    \item **allValues**: Vector of numerics for Mn values at each point.
#'  }
#' @export
#'
#' @examples
#' compute_Mn(electricity)
compute_Mn <- function(X, M = 10000, W = NULL, space = "BM", ...) {
  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }

  Zn <- (abs(.Zn(W,X, n)))^2
  return_value <- unname(unlist( 1 / M * rowSums(Zn) ))

  list(
    "value" = max(return_value),
    "location" = which.max(return_value),
    "allValues" = return_value
  )
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
.Zn <- function(W, X, n) {
  fhat_vals <- as.matrix(.fhat_all(X, W))

  sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])
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
.fhat_all <- function(X, W) {
  dot_col_cumsum(exp(1 / nrow(X) * complex(imaginary = 1) * t(X) %*% as.matrix(W)) / ncol(X))
}
