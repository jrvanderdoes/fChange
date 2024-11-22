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
#' @keywords internal
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
#' @keywords internal
#' @noRd
.fhat_all <- function(X, W) {
  dot_col_cumsum(exp(1 / nrow(X) * complex(imaginary = 1) * t(X) %*% as.matrix(W)) / ncol(X))
}
