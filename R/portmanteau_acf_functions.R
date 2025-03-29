#' Test Statistic - Welch Approximation for Q_WS_hyp_test
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
#' @param low_disc Boolean value specifying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#'         degrees of freedom (which approximates Q_h).
#'
#' @noRd
#' @keywords internal
Q_WS_quantile <- function(f_data, lag, alpha=0.05, M=NULL, low_disc=FALSE) {
  # Mean
  # mean_Q_h <- mean_hat_Q_h(f_data, lag)
  J <- NROW(f_data)
  f_data1 <- dfts(f_data)
  cov <- autocovariance(center(f_data1)$data^2,lags = lag, center=FALSE)

  mean_Q_h <- as.numeric(unlist(lapply(cov, function(x) sum(x))) / J^2)

  # Vars
  # var_Q_h <- variance_hat_Q_h(f_data, lag, M=M, low_disc=low_disc)
  var_Q_h <- sapply(lag,
                    function(lag,f_data, M, low_disc){
                      MCint_eta_approx_i_j(f_data, lag, lag, M=M, low_disc=low_disc)
                    },
                    f_data=f_data, M=M, low_disc=low_disc)

  # Get Welch Parameters
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * stats::qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag)
  p_val <- stats::pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, pvalue = p_val)
}


#' Test Statistic - Welch Approximation for iid Q_WS_hyp_test
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
#'
#' @noRd
#' @keywords internal
Q_WS_quantile_iid <- function(f_data, alpha=0.05) {
  J <- NROW(f_data)

  cov <- autocovariance(f_data,0)
  # Mean
  # mean_Q_h <- mean_hat_Q_h_iid(f_data)
  mean_Q_h <- (sum(diag(cov)) / J)^2

  # Var
  # var_Q_h <- variance_hat_Q_h_iid(f_data)
  var_Q_h <- 2 * ( sum(cov^2) / J^2 )^2

  # Welch Parameters
  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * mean_Q_h^2 / var_Q_h
  quantile <- beta * stats::qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag = 1)
  p_val <- stats::pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, pvalue = p_val)
}


#' Test Statistic - t_statistic_Q
#'
#' Computes the test statistic \code{Q_{T,h} = T*||y^hat_h||^2} for fixed h and for T
#'   inferred from the functional data f_data that is passed.
#'
#' @param f_data the functional data matrix with observed functions in columns
#' @param lag the fixed time lag used in the computation of the statistic
#'
#' @return scalar value of the statistic Q_{T,h} to test the hypothesis H_{0,h} : y_h(t,s) = 0.
#'
#' @noRd
#' @keywords internal
t_statistic_Q <- function(f_data, lag) {
  N <- NCOL(f_data)
  J <- NROW(f_data)

  # # Remove sapply
  # #gamma_hat <- autocov_approx_h(f_data, lag)
  # gamma_hat1 <- sapply(lag,
  #                      function(lag1,f_data){
  #                        # autocov_approx_h(f_data, lag1)
  #                        autocovariance(f_data,lag1)
  #                      },
  #                      f_data=f_data)
  #
  # #N * sum(gamma_hat^2) / (J^2)
  # N * colSums(gamma_hat1^2) / (J^2)

  if(length(lag)>1){
    val <- as.numeric(N * unlist(lapply(autocovariance(f_data,lag),function(x){sum(x^2)})) / J^2)
  }else{
    val <- N * sum(autocovariance(f_data,lag)^2) / J^2
  }

  val
}


#' Approximate MCint_eta_approx_i_j
#'
#' Computes an approximation of eta_i_j (defined under (15)) using the second
#'   Monte Carlo integration method "MCint" defined on page 8.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param i,j the indices i,j in 1:T that we are computing eta^hat_i_j for
#' @param M number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
#' @param low_disc boolean value specifying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar value of eta^_hat_i_j computed using the MCint method.
#'
#' @noRd
#' @keywords internal
MCint_eta_approx_i_j <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
  J <- NROW(f_data)
  N <- NCOL(f_data)

  if (is.null(M)) {
    M = floor((max(150 - N, 0) + max(100-J,0) + (J / sqrt(2))))
  }

  if (low_disc == TRUE) {
    stop('Low Discrepancy under work')
    # https://github.com/cran/fOptions/blob/master/R/LowDiscrepancy.R
    # rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
    # rand_samp_mat[which(rand_samp_mat == 0)] <- 1
  } else {
    rand_samp_mat <- matrix(nrow=M, ncol=4)
    for (k in 1:4) {
      rand_samp <- floor(J * stats::runif(M, 0, 1))
      rand_samp[which(rand_samp == 0)] <- 1
      rand_samp_mat[,k] <- rand_samp
    }
  }

  # eta_hat_i_j_sum <- 0
  # for (k in 1:M) {
  #   cov <- scalar_covariance_i_j(f_data, i, j, rand_samp_mat[k,])
  #   eta_hat_i_j_sum <- eta_hat_i_j_sum + cov^2
  # }
  ks <- (1+max(i,j)):N
  c_f_data <- center(f_data)

  tmp <- rowSums( c_f_data[rand_samp_mat[,1],ks-i,drop=FALSE] *
                    c_f_data[rand_samp_mat[,2],ks,drop=FALSE] *
                    c_f_data[rand_samp_mat[,3],ks-j,drop=FALSE] *
                    c_f_data[rand_samp_mat[,4],ks,drop=FALSE] ) / N
  eta_hat_i_j_sum <- sum(tmp^2)

  2/M * eta_hat_i_j_sum
}
