#' FACF - Welch Approximation
#'
#' Computes the 1-alpha quantile of the beta * chi-squared distribution with nu
#'   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#'   of the test statistic Q_h. This quantile is used to conduct an approximate size alpha test
#'   of the hypothesis H_0_h.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param lag specifies the lag used for the test statistic Q_h
#' @param alpha the significance level to be used in the hypothesis test
#' @param M optional argument specifying the sampling size in the related Monte Carlo method
#' @param low_disc boolean value specifiying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#'         degrees of freedom (which approximates Q_h).
Q_WS_quantile <- function(f_data, lag, alpha=0.05, M=NULL, low_disc=FALSE) {
  mean_Q_h <- mean_hat_Q_h(f_data, lag)
  var_Q_h <- variance_hat_Q_h(f_data, lag, M=M, low_disc=low_disc)

  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * stats::qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag)
  p_val <- stats::pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' FACF - Welch Approximation Mean Hat
#'
#' Computes the approximation of the mean defined in (15) which is used in the Welch-
#'   Satterthwaite approximation as mean of the chi-squared random variable approximating Q_h.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param lags specifies the lag use in the test statistic Q_h (lag = h in paper)
#'
#' @return scalar approximation of the mean of the test statistic Q_h.
mean_hat_Q_h <- function(f_data, lags) {
  J <- NROW(f_data)
  cov <- sapply(lags,
                function(lag,f_data){diagonal_covariance_i(f_data, lag)},
                f_data=f_data)
  colSums(cov) / J^2
}

#' Diagonal Covariance
#'
#' Computes the approximate diagonal covariance matrix of the functional
#'  data for lag windows defined by i.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i the index in 1:T that we are computing the covariance for
#'
#' @return A 2-dimensional array, encoding the covariance matrix of the functional data for lag
#' `windows defined by i.
diagonal_covariance_i <- function(f_data, i) {
  N <- NCOL(f_data)
  J <- NROW(f_data)
  c_f_data <- center(f_data)

  sum1 <- array(0, c(J, J))
  for (k in (1+i):N) {
    sum1 <- sum1 + ((c_f_data[,k-i])^2 %o% (c_f_data[,k])^2)
  }

  sum1 / N
  # c_f_data[,1:(N-i)]^2 %*% t(c_f_data[,(1+i):N]^2)/N
}


#' FACF - Welch Approximation Var Hat
#'
#' Computes the approximation of the variance defined in (15) which is used in
#'   the Welch- Satterthwaite approximation as variance of the chi-squared random variable
#'   approximating Q_h.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param lag specifies the lag use in the test statistic Q_h (lag = h in paper)
#' @param M optional argument specifying the sampling size in the related Monte Carlo method
#' @param low_disc boolean value specifiying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return Scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h <- function(f_data, lag, M=NULL, low_disc=FALSE) {
  sapply(lag,
         function(lag,f_data, M, low_disc){
           MCint_eta_approx_i_j(f_data, lag, lag, M=M, low_disc=low_disc)
         },
         f_data=f_data, M=M, low_disc=low_disc)
}


#' Monte Carlo Integration for eta_ij
#'
#' Computes an approximation of eta_i_j (defined under (15)) using the second
#'   Monte Carlo integration method "MCint" defined on page 8.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param i,j the indices in 1:T that we are computing eta^hat_i_j for
#' @param M number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
#' @param low_disc boolean value specifiying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar value of eta^_hat_i_j computed using the MCint method.
MCint_eta_approx_i_j <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
  J <- NROW(f_data)
  N <- NCOL(f_data)

  if (is.null(M)) {
    M = floor( max(150 - N, 0) + max(100-J, 0) + (J / sqrt(2)) )
  }
  if (low_disc == TRUE) {
    if (requireNamespace('fOptions', quietly = TRUE)) {
      rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
      rand_samp_mat[which(rand_samp_mat == 0)] <- 1
    } else {
      stop("Please install the 'fOptions' package for low discrepancy sampling.")
    }
  } else {
    rand_samp_mat <- matrix(nrow=M, ncol=4)
    for (k in 1:4) {
      rand_samp <- floor(J * stats::runif(M, 0, 1))
      rand_samp[which(rand_samp == 0)] <- 1
      rand_samp_mat[,k] <- rand_samp
    }
  }

  eta_hat_i_j_sum <- 0
  for (k in 1:M) {
    cov <- scalar_covariance_i_j(f_data, i, j, rand_samp_mat[k,])
    eta_hat_i_j_sum <- eta_hat_i_j_sum + cov^2
  }

  (2/M) * eta_hat_i_j_sum
}


#' Scalar Covariance for i,j
#'
#' Compute the approximate covariance at a point for lag windows defined by i,j
#'
#' \code{scalar_covariance_i_j} Computes the approximate covariance at a point of the functional data
#' for lag windows defined by i,j; a scalarized version of covariance_i_j that takes point estimates.
#' scalar_covariance_i_j returns the approximate covariance c^hat_i_j(t,s,u,v) evaluated at a
#'   given t,s,u,v in U_J X U_J X U_J X U_J (for use in MCint method).
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @param times A vector with 4 columns containing indices specifying which subset of f_data to consider (t,s,u,v)
#'
#' @return A numeric value; the covariance of the functional data at a point for lag
#'  windows defined by i,j (c^hat_i_j(t,s,u,v)).
scalar_covariance_i_j <- function(f_data, i, j, times) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  c_f_data <- center(f_data)

  sum1 <- 0
  for (k in (1+max(i,j)):N) {
    sum1 <- sum1 + c_f_data[times[1],k-i] * c_f_data[times[2],k] *
      c_f_data[times[3],k-j] * c_f_data[times[4],k]
  }
  # TODO:: Speed this up
  # ks <- (1+max(i,j)):N
  # sum2 <- c_f_data[times[1],ks-i] * c_f_data[times[2],ks] *
  #   c_f_data[times[3],ks-j] * c_f_data[times[4],ks]

  sum1 / N
}


#' FACF - Test Statistic
#'
#' Computes the test statistic \code{Q_{T,h} = T*||y^hat_h||^2} for fixed h and for T
#'   inferred from the functional data f_data that is passed.
#'
#' @param f_data the functional data matrix with observed functions in columns
#' @param lag the fixed time lag used in the computation of the statistic
#'
#' @return scalar value of the statistic Q_{T,h} to test the hypothesis H_{0,h} : y_h(t,s) = 0.
t_statistic_Q <- function(f_data, lag) {
  N <- NCOL(f_data)
  J <- NROW(f_data)

  #gamma_hat <- autocov_approx_h(f_data, lag)
  gamma_hat1 <- sapply(lag,
                       function(lag1,f_data){
                         autocov_approx_h(f_data, lag1)
                       },
                       f_data=f_data)

  #N * sum(gamma_hat^2) / (J^2)
  N * colSums(gamma_hat1^2) / (J^2)
}


#' Compute the approximate autocovariance at specified lag
#'
#' \code{autocov_approx_h} Computes the approximate autocovariance for a given lag h of the functional
#'  data for every (t,s) in U_J X U_J.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#'
#' @return A 2-dimensional array encoding the autocovariance matrix for a given lag h.
#' @export
autocov_approx_h <- function(f_data, lag) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)

  gamma_hat_sum <- 0
  for (i in 1:(N-lag)) {
    gamma_hat_sum <- gamma_hat_sum + c_f_data[,i] %o% c_f_data[,i+lag]
  }

  gamma_hat_sum / N
}


