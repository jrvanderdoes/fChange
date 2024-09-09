#' Block Bootstrapping
#'
#' Performs a block bootstrap on the functional data f_data with block size b.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param b the block size (of each block in each bootstrap sample)
#' @param B the number of bootstraps samples
#' @param moving boolean value specifying whether the block bootstrap should be moving or not. A moving black
#'  bootstrap samples individual functional observations and adds on the consequent block, rather than sampling
#'  blocks of the data.
#'
#' @return Returns a list of B elements, each element being a block bootstrap sample in the same format
#'  as the original functional data f_data.
#' @export
#'
#' @examples
#' res <- block_bootstrap(electricity,10)
block_bootstrap <- function(f_data, b, B = 300, moving = FALSE) {
  N <- NCOL(f_data)
  if (b > N) {
    stop("Please select a block size that is less than or equal to the sample size of
         the functional data. It is best to select a block size that evenly divides the
         sample size.")
  } else if (b < 1) {
    stop("The block size must be a positive integer.")
  } else if (B < 1) {
    stop("The number of bootstrap samples must be a positive integer.")
  }

  blocks <- list()
  M <- floor(N / b)
  for (s in 1:M) {
    blocks[[s]] <- (b*(s - 1) + 1):(b*s)
  }
  bootstrap_samples <- list()
  for (j in 1:B) {
    if (!moving) {
      samples <- sample(1:M, M, replace = TRUE)
      bootstrapped_data <- f_data[,blocks[[samples[1]]]]
      for (i in samples[-1]) {
        bootstrapped_data <- cbind(bootstrapped_data, f_data[,blocks[[samples[i]]]])
      }
    } else if (moving) {
      samples <- sample(1:(N - b), M, replace = TRUE)
      bootstrapped_data <- f_data[, samples[1]:(samples[1] + b)]
      for (i in 2:M) {
        bootstrapped_data <- cbind(bootstrapped_data, f_data[,samples[i]:(samples[i] + b)])
      }
    }
    bootstrap_samples[[j]] <- bootstrapped_data
  }

  bootstrap_samples
}

################################################################################

#' Compute the approximate autocovariance at specified lag
#'
#' \code{autocov_approx_h} Computes the approximate autocovariance for a given lag h of the functional
#' data
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#'
#' @return A 2-dimensional array encoding the autocovariance matrix for a given lag h.
autocov_approx_h <- function(f_data, lag) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)

  gamma_hat_sum <- 0
  for (i in 1:(N-lag)) {
    gamma_hat_sum <- gamma_hat_sum + c_f_data[,i] %o% c_f_data[,i+lag]
  }

  gamma_hat_sum / N
}


#' Compute the approximate covariance tensor for lag windows defined by i,j
#'
#' \code{covariance_i_j} Computes the approximate covariance tensor of the functional data for lag
#' windows defined by i,j.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#'
#' @return A 4-dimensional array, encoding the covariance tensor of the functional data for lag
#' windows defined by i,j.
covariance_i_j <- function(f_data, i, j) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)

  sum <- array(0, c(J, J, J, J))
  for (k in (1+max(i,j)):N) {
    sum <- sum + c_f_data[,k-i] %o% c_f_data[,k] %o% c_f_data[,k-j] %o% c_f_data[,k]
  }

  sum / N
}


# #' Compute the approximate covariance tensor for lag windows defined by i,j
# #'
# #' \code{covariance_i_j_vec} Computes the approximate covariance tensor of the functional data for lag
# #' windows defined by i,j; a vectorized version of covariance_i_j.
# #'
# #' @param f_data the functional data matrix with observed functions in the columns
# #' @param i,j the indices i,j in 1:T that we are computing the covariance for
# #' @return A 4-dimensional array, encoding the covariance tensor of the functional data for lag
# #' windows defined by i,j.
# covariance_i_j_vec <- function(f_data, i, j) {
#   N = NCOL(f_data)
#   J = NROW(f_data)
#
#   c_f_data <- center(f_data)
#   sum_parts <- as.list((1+max(i,j)):N)
#   sum_parts <- lapply(sum_parts,
#                       function(k) c_f_data[,k-i] %o% c_f_data[,k] %o%
#                         c_f_data[,k-j] %o% c_f_data[,k])
#
#   cov <- Reduce('+', sum_parts)
#   cov / N
# }


#' Compute the approximate diagonal covariance matrix for lag windows defined by i
#'
#' \code{diagonal_covariance_i} Computes the approximate diagonal covariance matrix of the functional
#' data for lag windows defined by i.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i the index in 1:T that we are computing the covariance for
#'
#' @return A 2-dimensional array, encoding the covariance matrix of the functional data for lag
#' windows defined by i.
diagonal_covariance_i <- function(f_data, i) {
  N = NCOL(f_data)
  J = NROW(f_data)
  c_f_data <- center(f_data)

  sum1 <- array(0, c(J, J))
  for (k in (1+i):N) {
    sum1 <- sum1 + ((c_f_data[,k-i])^2 %o% (c_f_data[,k])^2)
  }

  sum1 / N
}


#' Compute the approximate covariance at a point for lag windows defined by i,j
#'
#' \code{scalar_covariance_i_j} Computes the approximate covariance at a point of the functional data
#'  for lag windows defined by i,j; a scalarized version of covariance_i_j that takes point estimates.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' @param times A vector with 4 columns containing indices specifying which subset of f_data to consider
#'
#' @return A numeric value; the covariance of the functional data at a point for lag
#'  windows defined by i,j.
scalar_covariance_i_j <- function(f_data, i, j, times) {
  J <- NROW(f_data)
  N <- NCOL(f_data)
  c_f_data <- center(f_data)

  sum1 <- 0
  for (k in (1+max(i,j)):N) {
    sum1 <- sum1 + c_f_data[times[1],k-i] * c_f_data[times[2],k] * c_f_data[times[3],k-j]  *
      c_f_data[times[4],k]
  }

  sum1 / N
}


#' #' Compute the approximate covariance at a point for lag windows defined by i,j
#' #'
#' #' \code{scalar_covariance_i_j_vec} Computes the approximate covariance at a point of the functional data
#' #'  for lag windows defined by i,j; a vectorized version of \code{scalar_covariance_i_j}.
#' #'
#' #' @param f_data the functional data matrix with observed functions in the columns
#' #' @param i,j the indices i,j in 1:T that we are computing the covariance for
#' #' @param times A vector with 4 columns containing indices specifying which subset of f_data to consider
#' #'
#' #' @return A numeric value; the covariance of the functional data at a point for lag
#' #'  windows defined by i,j.
#' scalar_covariance_i_j_vec <- function(f_data, i, j, times) {
#'   J <- NROW(f_data)
#'   N <- NCOL(f_data)
#'   c_f_data <- center(f_data)
#'   sum_parts <- list((1+max(i,j)):N)
#'   sum_parts <- lapply(sum_parts,
#'                       function(k) c_f_data[times[1],k-i] * c_f_data[times[2],k] *
#'                         c_f_data[times[3],k-j]  * c_f_data[times[4],k])
#'   cov <- (1/N) * Reduce('+', sum_parts)
#'   cov
#' }


#' Compute part of the covariance under a strong white noise assumption
#'
#' \code{iid_covariance} A helper function used to compute one of the two independent sum terms in the
#' computation of the approximate covariance of the functional data under a strong white noise assumption.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @return A 2-dimensional matrix containing one of the two independent sums in the computation of the
#' covariance.
iid_covariance <- function(f_data) {
  N <- NCOL(f_data)
  c_f_data <- center(f_data)
  sum1 <- 0
  for (i in 1:N) {
    sum1 <- sum1 + c_f_data[,i] %o% c_f_data[,i]
  }
  sum1 / N
}


# #' Compute part of the covariance under a strong white noise assumption
# #'
# #' \code{iid_covariance_vec} A helper function used to compute one of the two independent sum terms in the
# #' computation of the approximate covariance of the functional data under a strong white noise assumption;
# #' a vectorized version of iid_covariance.
# #'
# #' @param f_data the functional data matrix with observed functions in the columns
# #' @return A 2-dimensional matrix containing one of the two independent sums in the computation of the
# #' covariance.
# iid_covariance_vec <- function(f_data) {
#   N <- NCOL(f_data)
#   c_f_data <- center(f_data)
#   sum_parts <- 1:N
#   sum_parts <- lapply(sum_parts,
#                       function(i) c_f_data[,i] %o% c_f_data[,i])
#   cov <- (1 / N) * Reduce('+', sum_parts)
#   cov
# }


#' List storage of diagonal covariances.
#'
#' \code{covariance_diag_store} Creates a list storage of approximate diagonal covariances computed
#' by the function \code{diagonal_covariance_i}.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param K the range of lags 1:K to use
#'
#' @return A list containing K 2-dimensional arrays containing the diagonal covariance matrices of the
#'  functional data, for lags h in the range 1:K.
covariance_diag_store <- function(f_data, K) {
  cov_i_store <- list()
  for (j in 1:K) {
    cov_i_store[[j]] <- diagonal_covariance_i(f_data, j)
  }
  cov_i_store
}


################################################################################

#' Test Statistic - Welch Approximation for V_WS_hyp_test
#'
#' Computes the 1-alpha quantile of the beta * chi-squared distribution with nu
#'   degrees of freedom, where beta and nu are obtained from a Welch-Satterthwaite approximation
#'   of the test statistic \code{V_K}. This quantile is used to conduct an approximate size alpha test
#'   of the hypothesis \code{H'_0_K}.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param K specifies the range of lags 1:K for the test statistic V_K
#' @param alpha the significance level to be used in the hypothesis test
#' @param M optional argument specifying the sampling size in the related Monte Carlo method
#' @param low_disc Boolean value specifying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return Scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#'           degrees of freedom (which approximates V_K)
V_WS_quantile <- function(f_data, K, alpha=0.05, M=NULL, low_disc=FALSE) {
  mean_V_K <- mean_hat_V_K(f_data, K)
  var_V_K <- variance_hat_V_K(f_data, K, M=M, low_disc=low_disc)

  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)

  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' Test Statistic - Welch Approximation for iid V_WS_hyp_test
#'
#' @inheritParams V_WS_quantile
#'
#' @return Scalar value of the 1-alpha quantile of the beta * chi-square distribution with nu
#'           degrees of freedom (which approximates V_K)
V_WS_quantile_iid <- function(f_data, K, alpha=0.05) {
  mean_V_K <- mean_hat_V_K_iid(f_data, K)
  var_V_K <- variance_hat_V_K_iid(f_data, K)

  beta <- var_V_K / (2 * mean_V_K)
  nu <- 2 * (mean_V_K^2) / var_V_K
  quantile <- beta * qchisq(1 - alpha, nu)

  statistic <- t_statistic_V(f_data, K)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


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
Q_WS_quantile <- function(f_data, lag, alpha=0.05, M=NULL, low_disc=FALSE) {
  mean_Q_h <- mean_hat_Q_h(f_data, lag)
  var_Q_h <- variance_hat_Q_h(f_data, lag, M=M, low_disc=low_disc)

  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
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
Q_WS_quantile_iid <- function(f_data, alpha=0.05) {
  mean_Q_h <- mean_hat_Q_h_iid(f_data)
  var_Q_h <- variance_hat_Q_h_iid(f_data)

  beta <- var_Q_h / (2 * mean_Q_h)
  nu <- 2 * (mean_Q_h^2) / var_Q_h
  quantile <- beta * qchisq(1 - alpha, nu)

  statistic <- t_statistic_Q(f_data, lag = 1)
  p_val <- pchisq(statistic / beta, nu, lower.tail = FALSE)

  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' Compute size alpha single-lag hypothesis test under weak or strong white noise assumption
#'
#' \code{Q_WS_hyp_test} Computes the size alpha test of a single lag hypothesis under a weak white noise
#' or strong white noise assumption using a Welch-Satterthwaite Approximation.
#'
#' @param f_data the functional data matrix with observed functions in the columns
#' @param lag the lag to use to compute the single lag test statistic
#' @param alpha the significance level to be used in the hypothesis test
#' @param iid boolean value, if given TRUE, the hypothesis test will use a strong-white noise assumption.
#' By default is FALSE, in which the hypothesis test will use a weak-white noise assumption.
#' @param M Number of samples to take when applying a Monte-Carlo approximation
#' @param low_disc Boolean value indicating whether or not to use low-discrepancy sampling in the Monte
#' Carlo method. Note, low-discrepancy sampling will yield deterministic results.
#' @param bootstrap boolean value, if given TRUE, the hypothesis test is done by approximating the
#' limiting distribution of the test statistic via a block bootstrap algorithm. FALSE by default
#' @param block_size the block size to be used in the block bootstrap method (in each bootstrap sample).
#' 10 by default.
#' @param straps the number of bootstrap samples to take; 300 by default
#' @param moving boolean value; determines whether or not the block bootstrap should be moving
#'
#' @return A list containing the p-value, the quantile, and a boolean value indicating whether or not the
#' hypothesis is rejected.
Q_WS_hyp_test <- function(f_data, lag, alpha=0.05, iid=FALSE,
                          M=NULL, low_disc=FALSE, bootstrap=FALSE,
                          block_size='adaptive', straps=300, moving = FALSE) {
  statistic <- t_statistic_Q(f_data, lag)

  if (bootstrap == TRUE) {
    if (block_size == 'adaptive') {
      block_size <- ceiling(adaptive_bandwidth(f_data, kernel = 'Bartlett'))
    }
    bootsraps <- list()
    bootstrap_samples <- block_bootsrap(f_data, block_size, B = straps, moving = moving)
    stats_distr <- lapply(bootstrap_samples, t_statistic_Q, lag=lag)
    statistic <- t_statistic_Q(f_data, lag=lag)
    quantile <- stats::quantile(as.numeric(stats_distr), 1 - alpha)
    p_value <- sum(statistic > stats_distr) / length(stats_distr)

    list(statistic = as.numeric(statistic), quantile = as.numeric(quantile),
         p_value = as.numeric(p_value), block_size = block_size)
  } else if (iid == FALSE) {
    results <- Q_WS_quantile(f_data, lag, alpha=alpha, M=M, low_disc=low_disc)
    statistic <- results$statistic
    quantile <- results$quantile
    p_val <- results$p_val
    reject <- statistic > quantile

    list(statistic = statistic, quantile = quantile, p_value = p_val)
  } else {
    results <- Q_WS_quantile_iid(f_data, alpha=alpha)
    statistic <- results$statistic
    quantile <- results$quantile
    p_val <- results$p_val
    reject <- statistic > quantile

    list(statistic= statistic, quantile = quantile, p_value = p_val)
  }
}

################################################################################

#' Compute Spectral Density Test
#'
#' Computes the spectral density operator based test statistic of the functional
#'   data f_data.
#'
#' @param f_data The functional data matrix with observed functions in columns
#' @param kernel The kernel function to use. The currently supported kernels are 'Bartlett' and 'Parzen'.
#                 The default kernel is 'Bartlett'.
#' @param bandwidth specifies the bandwidth to use. Currently admitted arguments are positive
#                    integers, 'static' which computes the bandwith p via p = n^(1/(2q+1)) where
#                    n is the sample size and q is the kernel order, or 'adaptive' which uses a
#                    bandwith selection method that is based on the functional data.
#'
#' @return scalar value of the spectral density based test statistic.
spectral_t_statistic <- function(f_data, kernel = 'Bartlett', bandwidth = 'adaptive') {
  J <- NROW(f_data)
  N <- NCOL(f_data)

  kernel_string <- kernel
  if (kernel == 'Bartlett') {
    kernel <- bartlett_kernel
    kernel_order <- 1
  } else if (kernel == 'Parzen') {
    kernel <- parzen_kernel
    kernel_order <- 2
    #} else if (kernel == 'Daniell') {
    #kernel <- daniell_kernel
    #kernel_order <- 2
  } else {
    stop("This kernel is not supported. Please see the documentation for supported kernel functions.")
  }

  if (bandwidth == 'static') {
    bandwidth <- N^(1 / (2 * kernel_order + 1))
  } else if (bandwidth == 'adaptive') {
    bandwidth <- max(2, adaptive_bandwidth(f_data, kernel_string))
  } else if (!is.numeric(bandwidth)) {
    stop("Please see the documentation for valid bandwith arguments.")
  }

  data_inner_prod <- crossprod(f_data) / J
  C_hat_HS_norm <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS_norm[j+1] <-
      N^(-2) * sum(data_inner_prod[(j+1):N, (j+1):N] *
                     data_inner_prod[1:(N-j), 1:(N-j)])
  }
  kernel_vals <- sapply(1:(N-1) / bandwidth, kernel)
  spectral_distance_Q_sq <- 2 * sum((kernel_vals^2) * C_hat_HS_norm[-1])
  C_n_k <- sum( (1 - 1:(N-1)/N) * kernel_vals^2 )
  D_n_k <- sum( (1 - 1:(N-2)/N) * (1 - 2:(N-1)/N) * kernel_vals[-(N-1)]^4 )
  sigma_squared_hat <- sum(diag(data_inner_prod)) / N

  t_stat_term <- sigma_squared_hat^(-2) * C_hat_HS_norm[1] * sqrt(2 * D_n_k) # for convenience
  untrans_num <- 2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq

  ### TODO: add case for when H is not R. This is denoted by const in original codebase
  beta <- 1 - (2/3) * sum(kernel_vals^2) * sum(kernel_vals^6) / (sum(kernel_vals^4)^2)
  t_stat <- ((2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq)^beta -
               (C_n_k^beta + 2^(-1) * beta * (beta - 1) * C_n_k^(beta-2) * t_stat_term^2)) / (beta * C_n_k^(beta-1) * t_stat_term)
  list(stat = t_stat, band = bandwidth)
}


#' Adaptive_bandwidth
#'
#' Computes the "optimal" bandwidth using a bandwidth selection method based on the
#'   spectral density operator which adapts to the functional data.
#'
#' @param f_data the functional data matrix with observed functions in columns
#' @param kernel the kernel function to use. The currently supported kernels are 'Bartlett' and 'Parzen'.
#'                 The default kernel is 'Bartlett'.
#'
#' @return a scalar value of the "optimal" data-adapted bandwidth.
adaptive_bandwidth <- function(f_data, kernel) {
  J <- NROW(f_data)
  N <- NCOL(f_data)

  if (kernel == 'Bartlett') {
    kernel <- bartlett_kernel
    order <- 1
    xi <- 1
    kern_int <- 2 / 3
  } else if (kernel == 'Parzen') {
    kernel <- parzen_kernel
    order <- 2
    xi <- 6
    kern_int <- 151 / 280
    #} else if (kernel == 'Daniell') {
    #kernel <- daniell_kernel
    #order <- 2
    #xi <- (pi^2) / 6
    #kern_int <- 1
  } else {
    stop('Please see the documentation for supported kernels.')
  }

  data_inner_prod <- crossprod(f_data) / (N * J)
  C_hat_HS <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS[j+1]<-sum(data_inner_prod[(j+1):N,(j+1):N] * data_inner_prod[1:(N-j),1:(N-j)])
  }

  initial_band_q <- 4 * N^(1 / (2*order +1))
  k_n_j_q <- kernel(1:(N-1) / initial_band_q)
  initial_band_0 <- initial_band_q / 4
  k_n_j_0 <- kernel(1:(N-1) / initial_band_0)

  Q_hat_sq <- 2 * sum(k_n_j_0^2 * C_hat_HS[-1])
  Q_hat_sq <- Q_hat_sq + sum(C_hat_HS[0])

  Term2 <- Q_hat_sq
  Term1 <- 2 * sum(k_n_j_q^2 * ((1:(N-1))^(2*order)) * C_hat_HS[-1])
  C_hat_TR <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_TR[j+1] <- sum(data_inner_prod[1:(N-j), (j+1):N])
  }
  trace <- 2 * sum(k_n_j_0^2 * C_hat_TR[-1])
  Term3 <- trace + sum(C_hat_TR[1])
  band_constant <- (2 * order * xi^2 * Term1 / (kern_int * (Term2 + Term3)))^(1/(2*order + 1))

  band_constant * N^(1 / (2*order + 1))
}


################################################################################

#' Test Statistic - t_statistic_Q
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


#' Test Statistic - t_statistic_V
#'
#' Computes the statistic \code{V_{T,K} = T*sum_h(||y^hat_h||^2)} or h in 1:K and for T
#'   inferred from the functional data f_data that is passed to the function.
#'
#' @param f_data the functional data with functions in columns
#' @param K the max value in the range of time lags (1:K) used
#'
#' @return scalar value of the statistic V_{T,K} to test the hypothesis
#'         H'_{0,K} : for all h in 1:K y_h(t,s) = 0.
t_statistic_V <- function(f_data, K) {
  V_T_K <- 0
  for (h in 1:K) {
    V_T_K <- V_T_K + t_statistic_Q(f_data, h)
  }

  V_T_K
}

################################################################################

#' Approximate eta_i_j
#'
#' a non-stochaistic approximation of eta_i_j (see (15)) using a Riemann sum.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param i,j the indices i,j in 1:T that we are computing eta^hat_i_j for
#'
#' @return scalar value of eta^hat_i_j computed using a simple Riemann sum.
true_eta_approx_i_j <- function(f_data, i, j) {
  J <- NROW(f_data)
  cov_tensor <- covariance_i_j(f_data, i, j)
  2 * sum(cov_tensor^2) / J^4
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
MCint_eta_approx_i_j <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
  J <- NROW(f_data)
  N <- NCOL(f_data)

  if (is.null(M)) {
    M = floor((max(150 - N, 0) + max(100-J,0) + (J / sqrt(2))))
  }

  if (low_disc == TRUE) {
    if (requireNamespace('fOptions', quietly = TRUE)) {
      rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
      rand_samp_mat[which(rand_samp_mat == 0)] = 1
    } else {
      stop("Please install the 'fOptions' package for low discrepancy sampling.")
    }
  } else {
    rand_samp_mat <- matrix(nrow=M, ncol=4)
    for (k in 1:4) {
      rand_samp <- floor(J * runif(M, 0, 1))
      rand_samp[which(rand_samp == 0)] = 1
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

# # MCint_eta_approx_i_j_vec is a vectorized version of MCint_eta_approx_i_j.
# # Input: f_data = the functional data matrix with functions in columns
# #        i,j = the indices i,j in 1:T that we are computing eta^hat_i_j for
# #        M = number of vectors (v1, v2, v3, v4) to sample uniformly from U_J X U_J X U_J X U_J
# #        low_disc = boolean value specifiying whether or not to use low-discrepancy sampling
# #                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
# # Output: scalar value of eta^_hat_i_j computed using the MCint method.
# MCint_eta_approx_i_j_vec <- function(f_data, i, j, M=NULL, low_disc=FALSE) {
#   J <- NROW(f_data)
#   N <- NCOL(f_data)
#   M = floor((max(150 - N, 0) + max(100-J,0) + (J / sqrt(2))))
#   if (low_disc == TRUE) {
#     if (requireNamespace('fOptions', quietly = TRUE)) {
#       rand_samp_mat <- apply(J * fOptions::runif.sobol(M, 4, scrambling = 3), 2, floor)
#       rand_samp_mat[which(rand_samp_mat == 0)] = 1
#     } else {
#       stop("Please install the 'fOptions' package for low discrepancy sampling.")
#     }
#   } else {
#     rand_samp_mat <- matrix(nrow=M, ncol=4)
#     for (k in 1:4) {
#       rand_samp <- floor(J * runif(M, 0, 1))
#       rand_samp[which(rand_samp == 0)] = 1
#       rand_samp_mat[,k] <- rand_samp
#     }
#   }
#   eta_parts <- as.list(1:M)
#   eta_parts <- lapply(eta_parts, function(k) scalar_covariance_i_j(f_data, i, j,
#                                                                    rand_samp_mat[k,]) ^ 2)
#   eta_hat_i_j <- (2 / M) * Reduce('+', eta_parts)
#   eta_hat_i_j
# }


#' Mean Approximation - mean_hat_V_K
#'
#' Computes the approximation of the mean defined in (15) which is used in the Welch-
#'   Satterthwaite approximation as mean of the chi-squared random variable approximating V_K.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param K specifies the range of lags 1:K for for the test statistic V_K
#'
#' @return scalar approximation of the mean of the test statistic V_K.
mean_hat_V_K <- function(f_data, K) {
  J <- NROW(f_data)
  sum1 <- 0
  store <- covariance_diag_store(f_data, K)
  for (i in 1:K) {
    sum1 <- sum1 + sum(store[[i]])
  }

  sum1 / J^2
}


#' Mean Approximation - iid mean_hat_V_K
#'
#' Computes the approximation of the mean defined in (15) which is used in the
#'   Welch-Satterthwaite approximation under the assumption that the functional data follows a
#'   strong white noise.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param K specifies the range of lags 1:K for for the test statistic V_K
#'
#' @return scalar approximation of the mean of the test statistic V_K under a strong white noise
#'         assumption.
mean_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov <- iid_covariance(f_data)

  K * ( (sum(diag(cov)) / J)^2)
}


#' Mean Approximation - mean_hat_Q_h
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


#' Mean Approximation - iid mean_hat_Q_h
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


#' Variance Approximation - variance_hat_V_K
#'
#' Computes the approximation of the variance defined in (15) which is used in
#'   the Welch- Satterthwaite approximation as the variance of the chi-squared random variable
#'   approximating V_K.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param K specifies the range of lags 1:K for the test statistic V_K
#' @param M optional argument specifying the sampling size in the related Monte Carlo method
#' @param low_disc boolean value specifiying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar approximation of the variance of the test statistic V_K.
variance_hat_V_K <- function(f_data, K, M=NULL, low_disc=FALSE) {
  N <- NCOL(f_data)
  sum1 <- 0
  for (i in 1:K) {
    sum1 <- sum1 + MCint_eta_approx_i_j(f_data, i, i, M=M, low_disc=low_disc)
  }

  bandwidth <- ceiling(0.25 * (N ^ (1/3)))
  if (K > 1) {
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        if (abs(i-j) > bandwidth) { # empirically, past a lag of 15, error is less than 1%
          next
        }
        sum1 <- sum1 + (2 * MCint_eta_approx_i_j(f_data, i, j, M=M, low_disc=low_disc))
      }
    }
  }

  sum1
}


#' Variance Approximation - iid variance_hat_V_K
#'
#' Computes the approximation of the variance defined in (15) which is used
#'   in the Welch-Satterthwaite approximation under the assumption that the functional data
#'   follows a strong white noise.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param K specifies the range of lags 1:K for the test statistic V_K
#'
#' @return scalar approximation of the variance of the test statistic V_K
variance_hat_V_K_iid <- function(f_data, K) {
  J <- NROW(f_data)
  cov_iid <- iid_covariance(f_data)

  K * 2 * ( sum(cov_iid^2) / (J^2) )^2
}


#' Variance Approximation - variance_hat_Q_h
#'
#' Computes the approximation of the variance defined in (15) which is used in
#'   the Welch-Satterthwaite approximation as variance of the chi-squared random variable
#'   approximating Q_h.
#'
#' @param f_data the functional data matrix with functions in columns
#' @param lag specifies the lag use in the test statistic Q_h (lag = h in paper)
#' @param M optional argument specifying the sampling size in the related Monte Carlo method
#' @param low_disc boolean value specifiying whether or not to use low-discrepancy sampling
#'                   for the Monte-Carlo method (only Sobol Sampling is currently supported)
#'
#' @return scalar approximation of the variance of the test statistic Q_h
variance_hat_Q_h <- function(f_data, lag, M=NULL, low_disc=FALSE) {
  sapply(lag,
         function(lag,f_data, M, low_disc){
           MCint_eta_approx_i_j(f_data, lag, lag, M=M, low_disc=low_disc)
         },
         f_data=f_data, M=M, low_disc=low_disc)
}


#' Variance Approximation - iid variance_hat_Q_h
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

