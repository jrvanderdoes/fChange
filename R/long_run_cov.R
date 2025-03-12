#' Estimate Long-run Covariance Kernel
#'
#' This (internal) function estimates the long-run covariance kernel. That is,
#'     \eqn{C_{\epsilon}(t,t') = \sum_{l=-\inf}^{\inf} \text{Cov}(\epsilon_0(t),
#'     \epsilon_l(t'))} with error sequence \eqn{(\epsilon_i : i \in \mathbb{Z})}.
#'
#' This is an internal function, see use in .change_mean.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param  h The window parameter parameter for the estimation of the long run
#'  covariance kernel. The default value is \code{h=0}, i.e., it assumes iid
#'  data. Note there exists an internal check such that \eqn{h=min(h,ncol(X)-1)}.
#' @param K (Optional) Function indicating the Kernel to use if h>0
#'
#' @return Data.frame of numerics with dim of ncol(data) x ncol(data), that is
#'     symmetric.
#' @export
#'
#' @examples
#' result <- long_run_covariance(electricity,2)
long_run_covariance <- function(X, h=0, K=bartlett_kernel){
  X <- center(dfts(X))

  N <- ncol(X$data)
  res <- nrow(X$data)
  Ceps <- matrix(NA, nrow = res, ncol = res)

  h <- min(h, N-1)

  for (k in 1:res) {
    for (r in k:res) {
      # Multiple all observations taken at the same point in time across FDs
      s <- X$data[k, ] %*% X$data[r, ]
      if (h > 0) {
        # Take care of lagged values if important
        for (i in 1:h) {
          # if( (N - i)>=1 && (i + 1)<=N){}
          a <- as.numeric(X$data[k, 1:(N - i)]) %*% as.numeric(X$data[r, (i + 1):N])
          a <- a + as.numeric(X$data[r, 1:(N - i)]) %*% as.numeric(X$data[k, (i + 1):N])
          s <- s + K(i/h) * a
        }
      }
      Ceps[k, r] <- Ceps[r, k] <- s
    }
  }

  Ceps / N
}
