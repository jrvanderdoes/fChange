#' Estimate Long-run Covariance Kernel
#'
#' This (internal) function estimates the long-run covariance kernel. That is,
#'     \eqn{C_{\epsilon}(t,t') = \sum_{l=-\inf}^{\inf} \text{Cov}(\epsilon_0(t),
#'     \epsilon_l(t'))} with error sequence \eqn{(\epsilon_i : i \in \mathbb{Z})}.
#'
#' This is an internal function, see use in .change_mean.
#'
#' @param data funts object or Numeric data.frame with evaled points on rows and fd objects in columns
#' @param  h The window parameter parameter for the estimation of the long run covariance kernel. The default
#'  value is \code{h=0}, i.e., it assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#'
#' @return Data.frame of numerics with dim of ncol(data) x ncol(data), that is
#'     symmetric.
#'
#' @noRd
#' @keywords internal
.long_run_cov <- function(data, h, K){
  data <- center(funts(data))

  N <- ncol(data$data)
  res <- nrow(data$data)
  Ceps <- matrix(NA, nrow = res, ncol = res)

  for (k in 1:res) {
    for (r in k:res) {
      # Multiple all observations taken at the same point in time across FDs
      s <- data$data[k, ] %*% data$data[r, ]
      if (h > 0) {
        # Take care of lagged values if important
        for (i in 1:h) {
          a <- as.numeric(data$data[k, 1:(N - i)]) %*% as.numeric(data$data[r, (i + 1):N])
          a <- a + as.numeric(data$data[r, 1:(N - i)]) %*% as.numeric(data$data[k, (i + 1):N])
          s <- s + K(i/h) * a
        }
      }
      Ceps[k, r] <- Ceps[r, k] <- s
    }
  }

  Ceps / N
}
