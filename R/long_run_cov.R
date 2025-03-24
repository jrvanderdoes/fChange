#' Estimate Long-run Covariance Kernel
#'
#' Estimate the long-run covariance kernel for functional data. That is, solve
#'  \eqn{
#'    C_{\epsilon}(t,t') = \sum_{l=-\inf}^{\inf} \text{Cov}(\epsilon_0(t),
#'     \epsilon_l(t'))
#'  }
#'  with sequence \eqn{(\epsilon_i : i \in \mathbb{Z})} defined as the centered
#'  data (can center based on changes if given).
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param  h The window parameter parameter for the estimation of the long run
#'  covariance kernel. The default value is \code{h=2*ncol(X)^(1/5)}.
#'  Note there exists an internal check such that \eqn{h=min(h,ncol(X)-1)} when
#'  alternative options are given.
#' @param K Function indicating the kernel to use if \eqn{h>0}.
#' @param changes Vector of numeric change point locations. Can be NULL.
#'
#' @return Symmetric data.frame of numerics with dim of ncol(data) x ncol(data).
#' @export
#'
#' @examples
#' result <- long_run_covariance(electricity,2)
long_run_covariance <- function(X, h=2*ncol(X)^(1/5), K=bartlett_kernel,
                                changes=NULL){
  X <- center(dfts(X), changes = changes)

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
