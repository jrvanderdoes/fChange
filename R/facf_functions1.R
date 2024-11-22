#' FACF - Test Statistic H_0
#'
#' Computes the size alpha test of the hypothesis H_0_h using the WS
#'   Approximation under the assumption that the data follows a strong white noise.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param alpha the significance level to be used in the hypothesis test
#'
#' @return scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#'         degrees of freedom (which approximates Q_h) (computed under a strong white noise
#'         assumption).
Q_WS_quantile_iid <- function(f_data, alpha=0.05) {
  mean_Q_h <- mean_hat_Q_h_iid(f_data)
  var_Q_h <- variance_hat_Q_h_iid(f_data)

  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * stats::qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag = 1)
  p_val <- stats::pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' Test Statistic - Mean Hat Q_h_iid
#'
#' Computes the approximation of the mean defined in (15) which is used in the
#'   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#'   strong white noise.
#'
#' @param f_data the functional data matrix with functions in columns
#'
#' @return scalar approximation of the mean of the test statistic Q_h under a strong white noise
#'         assumption.
mean_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)

  (sum(diag(cov)) / J)^2
}


#' Partial Covariance under iid
#'
#' Compute part of the covariance under a strong white noise assumption
#'
#' \code{iid_covariance} A helper function used to compute one of the two independent sum terms in the
#' `computation of the approximate covariance of the functional data under a strong white noise assumption.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#'
#' @return A 2-dimensional matrix containing one of the two independent sums in the computation of the
#'  covariance.
iid_covariance <- function(f_data) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  # sum1 <- 0
  # for (i in 1:N) {
  #   sum1 <- sum1 + c_f_data[,i] %o% c_f_data[,i]
  # }
  # sum1 / N
  c_f_data[,1:N] %*% t(c_f_data[,1:N]) / N
}


#' Test Statistic - Variance Q_h_iid
#'
#' Computes the approximation of the variance defined in (15) which is used
#'   in the Welch- Satterthwaite approximation under the assumption that the functional data
#'   follows a strong white noise.
#'
#' @param f_data the functional data matrix with functions in columns
#'
#' @return scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h_iid <- function(f_data) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)

  2 * ( sum(cov_iid^2) / (J^2) )^2
}
