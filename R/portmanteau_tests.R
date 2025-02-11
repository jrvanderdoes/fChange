#' Compute Functional Hypothesis Tests
#'
#' Computes a variety of functional portmanteau hypothesis tests. All hypothesis tests in this
#'  package are accessible through this function.
#'
#' @details The "single-lag" portmanteau test is based on the sample autocovariance function computed from the
#' functional data. This test assesses the significance of lagged autocovariance operators at a single, user-specified
#' lag h. More specifically, it tests the null hypothesis that the lag-h autocovariance operator is equal to 0.
#' This test is designed for stationary functional time-series, and is valid under conditional heteroscedasticity
#' conditions.
#' The required parameter for this test are 'lag', which determines the lag at which the test is evaluated. If this
#' parameter is left blank, it will take a default of 1.
#' The optional parameters for this test are 'iid', 'M', 'low_disc', 'bootstrap', 'block_size', 'straps', 'moving',
#' and 'alpha'.
#'
#' The "multi-lag" portmanteau test is also based on the sample autocovariance function computed from the functional
#' data. This test assesses the cumulative significance of lagged autocovariance operators, up to a user-selected
#' maximum lag lag. More specifically, it tests the null hypothesis that the first lag lag-h autocovariance operators
#' (h going from 1 to lag) is equal to 0. This test is designed for stationary functional time-series, and is valid
#' under conditional heteroscedasticity conditions.
#' The required parameter for this test is 'lag', which determines the maximum lag at which the test is evaluated.
#' If this parameter is left blank, it will take a default of 20.
#' The optional parameters for this test are 'iid', 'M', 'low_disc', 'bootstrap', 'block_size', 'straps', 'moving',
#' and 'alpha'.
#'
#' The "spectral" portmanteau test is based on the spectral density operator. It essentially measures the proximity of a
#' functional time series to a white noise - the constant spectral density operator of an uncorrelated series.
#' Unlike the "single-lag" and "multi-lag" tests, this test is not for general white noise series, and may not hold
#' under functional conditionally heteroscedastic assumptions.
#' The optional parameters for this test are 'kernel', 'bandwidth', and 'alpha'.
#'
#' The "independence" portmanteau test is a test of independence and identical distribution based on a dimensionality
#' reduction by projecting the data onto the most important functional principal components. It is based on the
#' resulting lagged cross-variances. This test is not for general white noise series, and may not hold under
#' functional conditionally heteroscedastic assumptions.
#' The required parameters for this test are 'lag' and 'components'. The 'lag' parameter determines the maximum lag at
#' which the test is evaluated. The 'components' parameter determines the number of the most important principal
#' components to use (importance is determined by the proportion of the variance that is explained by the
#' individual principal component.)
#'
#' The "imhof" portmanteau test is an analogue of the "single-lag" test. While the "single-lag" test computes the
#' limiting distribution of the test statistic via a Welch-Satterthwaite approximation, the "imhof" test directly
#' computes the coefficients of the quadratic form in Normal variables which the test statistic converges too as
#' the sample size goes to infinity. We warn the user that this test is extremely computationally expensive, and
#' is only recommended for small datasets as a means of cross-verification against the single-lag test.
#' The required parameter for this test is 'lag', which determines the lag at which the test is evaluated.
#' The "imhof" test requires the "tensorA" and "CompQuadForm" packages. Note also that the imhof test does not
#' return a statistic, and thus returns a list with only 2 elements if suppress_raw_output = FALSE.
#'
#' @param data A funts object or functional data matrix with observed functions in the columns.
#' @param test A String specifying the hypothesis test. Currently available tests are referred to by their
#'  string handles: 'variety', "single-lag", "multi-lag", "spectral", "independence", and "imhof". Please see the Details
#'  section of the documentation, or the vignette, for a short overview of the available tests. For a more
#'  complete treatment of these hypothesis tests, please consult the references.
#' @param lag A positive integer value. Only used for the "single-lag", "multi-lag", "independence", and "imhof" tests.
#'  This parameter specifies the single lag, or maximum lag, to be used by the specified test.
#' @param M Only used for the "single-lag" and "multi-lag" tests. A positive Integer. Determines the number of
#'  Monte-Carlo simulations employed in the Welch-Satterthwaite approximation of the limiting distribution of the
#'  test statistic.
#' @param method String indicating the method, options include:
#'  \itemize{
#'    \item **iid**:Only used for the "single-lag" and "multi-lag" tests. The
#'      hypothesis test will use a strong-white noise assumption (instead of a
#'      weak-white noise assumption).
#'    \item **bootstrap**: Only used for the "single-lag" test. The hypothesis
#'      test is evaluated by approximating the limiting distribution of the test
#'      statistic via a block bootstrapping process.
#'  }
#' @param kernel Only used for the "spectral" test. A String, 'Bartlett' by default. Specifies the kernel to be
#'  used in the "spectral" test. Currently supported kernels are the 'Bartlett' and 'Parzen' kernels.
#' @param bandwidth Only used for the "spectral" test. Either a String or a positive Integer value, 'adaptive' by
#'  default. Determines the bandwidth (or lag-window) to be used for the test. Given the string handle 'adaptive',
#'  the bandwidth is computed via a bandwidth selection method which aims to minimize the integrated normed
#'  error of the spectral density operator. If the given string handle is 'static', the bandwidth is computed
#'  to be n^(1/(2q + 1)), where n is the sample size and q is the kernel order. If a positive integer is
#'  given, that will be the bandwidth that is used.
#' @param components Only used for the "independence" test. A positive Integer value. Determines the number of
#'  functional principal components to use (ranked by their importance).
#' @param block_size Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A positive Integer
#'  value, with the default value being computed via the adaptive bandwidth selection method in the "spectral" test.
#'  Determines the block size (of each block in each bootstrap sample) if the test is being bootstrapped.
#' @param moving Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A Boolean value, FALSE
#'  by default If given TRUE, the performed block bootstrap will be moving rather than stationary.
#' @param B Only used for the "single-lag" test in the case when 'bootstrap' = TRUE. A positive Integer with
#'  a default value of 300. Determines the number of bootstrap samples to take if the test is being bootstrapped.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#'  hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#'  1-alpha quantile of the limiting distribution of the specified test's test statistic.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#'  limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#'  containing a short description of the test, the p-value, and additional information about the test if
#'  suppress_print_output = FALSE. If 'complete-test' = TRUE, will return a 1-column table instead containing
#'  the p-values for a variety of tests, which are given short descriptions in the index of the table.
#' @export
#'
#' @references Kim, M., Kokoszka P., & Rice G. (2023) White noise testing for functional
#'  time series. Statist. Surv., 17, 119-168, DOI: 10.1214/23-SS143
#'
#' @references Characiejus V., & Rice G. (2019). A general white noise test based on kernel lag-window estimates of the
#' spectral density operator. Econometrics and Statistics, submitted.
#'
#' @references Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @references Zhang X. (2016). White noise testing and model diagnostic checking for functional time series.
#' Journal of Econometrics, 194, 76-95.
#'
#' @references Gabrys R., & Kokoszka P. (2007). Portmanteau Test of Independence for Functional Observations.
#' Journal of the American Statistical Association, 102:480, 1338-1348, DOI: 10.1198/016214507000001111.
#'
#' @references Chen W.W. & Deo R.S. (2004). Power transformations to induce normality and their applications.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66, 117-130.
#'
#' @examples
#' b <- generate_brownian_motion(250)
#' res0 <- portmanteau_tests(b, test = 'single-lag', lag = 10)
#' res1 <- portmanteau_tests(b, test = 'multi-lag', lag = 10, alpha = 0.01)
#' res2 <- portmanteau_tests(b, test = 'spectral', kernel = 'Bartlett',
#'                           bandwidth = 'static', alpha = 0.05)
#' res3 <- portmanteau_tests(b, test = 'spectral', alpha = 0.1,
#'                           kernel = 'Parzen', bandwidth = 'adaptive')
#' res4 <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
#' res5 <- portmanteau_tests(b, test = 'single-lag', lag = 1, M = 250)
portmanteau_tests <- function(data, test = c('variety', 'single-lag', 'multi-lag',
                                      'spectral', 'independence', 'imhof'),
                       lag=NULL, M=NULL,
                       method = c('iid','bootstrap'),
                       kernel = "Bartlett", bandwidth = "adaptive",
                       components = 3, block_size = "adaptive", moving=FALSE,
                       B = 500, alpha=0.05) {
  data <- funts(data)
  poss_tests <- c('variety', 'single-lag', 'multi-lag', 'spectral',
                    'independence', 'imhof')
  test <- .verify_input(test,poss_tests)
  method <- .verify_input(method, c('iid','bootstrap'))

  # Alpha
  if (alpha < 0 | alpha > 1) {
    stop("Invalid arguments, the significance level alpha must be between 0 and 1.")
  }

  if (!is.null(lag)) {
    if (!all.equal(lag, as.integer(lag)) | lag <= 0) {
      stop("Invalid arguments, lag must be a positive integer for the single-lag and multi-lag tests.")
    }
  }
  if (!is.null(M)) {
    if (!all.equal(M, as.integer(M)) | M < 0) {
      stop("Invalid arguments, M must be a positive integer or NULL.")
    }
  }

  ## RUN TESTS
  if (test == 'variety') {
    m <- as.table(matrix(0, ncol = 1, 10))
    colnames(m) <- c('p_value')
    rownames(m) <- c('single-lag, lag = 1', 'single-lag, lag = 2',
                     'single-lag, lag = 3', 'multi-lag, lag = 5',
                     'multi-lag, lag = 10', 'multi-lag, lag = 20',
                     'spectral, static bandwidth',
                     'spectral, adaptive bandwidth',
                     'independence, 3 components, lag = 3',
                     'independence, 16 components, lag = 10')

    m[1] <- portmanteau_tests(data, test = 'single-lag', lag = 1, alpha=alpha)$p_value
    m[2] <- portmanteau_tests(data, test = 'single-lag', lag = 2, alpha=alpha)$p_value
    m[3] <- portmanteau_tests(data, test = 'single-lag', lag = 3, alpha=alpha)$p_value

    m[4] <- portmanteau_tests(data, test = 'multi-lag', lag = 5, alpha=alpha)$p_value
    m[5] <- portmanteau_tests(data, test = 'multi-lag', lag = 10, alpha=alpha)$p_value
    m[6] <- portmanteau_tests(data, test = 'multi-lag', lag = 20, alpha=alpha)$p_value

    m[7] <- portmanteau_tests(data, test = 'spectral', method=method, B=B,
                       bandwidth = 'static')$p_value
    m[8] <- portmanteau_tests(data, test = 'spectral', method=method, B=B,
                       bandwidth = 'adaptive')$p_value

    m[9] <- portmanteau_tests(data, test = 'independence', method=method, B=B,
                       components = 3, lag = 3)$p_value

    m[10] <- portmanteau_tests(data, test = 'independence', method=method, B=B,
                        components = 16, lag = 10)$p_value
    results <- m
  } else if (test == 'multi-lag') {
    # Lag Tests
    if (is.null(lag)) {
      warning("No maximal lag given. Default to lag = 20.")
      lag = 20
    }

    results <- .multi_lag_test(data, lag, M=M, method=method, alpha=alpha)
  } else if (test == 'single-lag') {
    if (is.null(lag)) {
      warning("No maximal lag given. Default to lag = 1.")
      lag = 1
    }
    results <- .single_lag_test(data, lag, alpha=alpha,
                    M=M, method=method,
                    block_size=block_size, B=B,
                    moving = moving)
  } else if (test == 'spectral') {
    results <- .spectral_test(data, kernel = kernel, bandwidth = bandwidth, alpha = alpha)
  } else if (test == 'independence') {
    results <- .independence_test(data, components = components, lag = lag, alpha = alpha)
  } else if (test == 'imhof') {
    input <- readline("The imhof test is computationally expensive. \n
                      Press [enter] if you would like to continue.")
    if (input != '') {
      stop("User cancelled the test.")
    }
    results <- .imhof_test(data, lag)
  }

  results
}


#' Single-Lag Hypothesis Test
#'
#' Computes the single-lag hypothesis test at a single user-specified lag.
#'
#' @details The "single-lag" portmanteau test is based on the sample autocovariance function computed from the
#'  functional data. This test assesses the significance of lagged autocovariance operators at a single,
#'  user-specified lag h. More specifically, it tests the null hypothesis that the lag-h autocovariance
#'  operator is equal to 0. This test is designed for stationary functional time-series, and is valid under
#'  conditional heteroscedasticity conditions.
#'
#' @param data funts object or functional data matrix with observed functions in the columns
#' @param lag Positive integer value. The lag to use to compute the single lag test statistic.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#'  hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#'  1-alpha quantile of the limiting distribution of the specified test's test statistic.
#' @param method string indicating type of test.
#'  \itemize{
#'    \item **iid**: The hypothesis test will use a strong-white
#'      noise assumption (instead of a weak-white noise assumption).
#'    \item **lowdiscrepancy**: The hypothesis test will usea
#'      low-discrepancy sampling in the Monte-Carlo method. Note,
#'      low-discrepancy sampling will yield deterministic results.
#'    \item **bootstrap**: The hypothesis test is done by approximating the
#'      limiting distribution of the test statistic via a block bootstrap
#'      process.
#'  }
#' @param M Positive integer value. Number of Monte-Carlo simulations for the Welch-Satterthwaite approximation.
#' @param block_size A positive Integer value, with the default value being computed via the adaptive
#'  bandwidth selection method in the "spectral" test. Determines the block size (of each block in each
#'  bootstrap sample) if the test is being bootstrapped.
#' @param B A positive Integer, with a default value of 300. Determines the number of bootstrap samples
#'  to take if the test is being bootstrapped. Only used if 'bootstrap' == TRUE.
#' @param moving A Boolean value, FALSE by default If given TRUE, the performed block bootstrap will be moving
#'  rather than stationary.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#'  limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#'  containing a short description of the test, the p-value, and additional information about the test if
#'  suppress_print_output = FALSE.
#'
#' @noRd
#' @keywords internal
#'
#' @references Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#' under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' # f <- generate_brownian_motion(100)
#' # .single_lag_test(f, lag = 1)
#' # .single_lag_test(f, lag = 2, M=100)
.single_lag_test <- function(
    data, lag=1, alpha=0.05, method = c('iid', 'lowdiscrepancy', 'bootstrap'),
    M=NULL, block_size='adaptive', B=300, moving = FALSE){
  data <- funts(data)

  poss_methods <- c('iid', 'lowdiscrepancy', 'bootstrap')
  method <- .verify_input(method, poss_methods)

  if(method=='bootstrap') {
    if (block_size == 'adaptive') {
      block_size <- ceiling(adaptive_bandwidth(data$data, kernel = 'Bartlett'))
    }

    # bootstrap_samples <- block_bootstrap(data$data, block_size, B = B, moving = moving)
    stats_distr <- .bootstrap(data, blocksize = block_size,
                                    M = B, type =  ifelse(moving,'overlapping', 'separate'),
                                    replace = TRUE,lag=lag, fn=t_statistic_Q)
    # stats_distr <- lapply(bootstrap_samples, t_statistic_Q, lag=lag)
    quantile <- stats::quantile(as.numeric(stats_distr), 1 - alpha)
    statistic <- t_statistic_Q(data$data, lag)
    p_value <- sum(statistic <= stats_distr) / length(stats_distr)

    results <- list(statistic = as.numeric(statistic),
                    quantile = as.numeric(quantile),
                    p_value = as.numeric(p_value), block_size = block_size)
  } else if(method=='lowdiscrepancy') {
    results <- Q_WS_quantile(data$data, lag, alpha=alpha, M=M, low_disc=TRUE)
  } else if(method=='iid') {
    results <- Q_WS_quantile_iid(data$data, alpha=alpha)
  }

  results
}


#' Multi-Lag Hypothesis Test
#'
#' Computes the multi-lag hypothesis test over a range of user-specified lags.
#'
#' @details The "multi-lag" portmanteau test is also based on the sample autocovariance function computed from the
#'  functional data. This test assesses the cumulative significance of lagged autocovariance operators, up to a
#'  user-selected maximum lag lag. More specifically, it tests the null hypothesis that the first lag lag-h autocovariance
#'  operators (h going from 1 to lag) is equal to 0. This test is designed for stationary functional time-series, and
#'  is valid under conditional heteroscedasticity conditions.
#'
#' @param data funts object or functional data matrix with observed functions in the columns
#' @param lag Positive integer value. The lag to use to compute the single lag test statistic
#' @param M Positive integer value. Number of Monte-Carlo simulation for Welch-Satterthwaite approximation.
#' @param method Method to use:
#'  \itemize{
#'    \item **iid**: The hypothesis test will use a strong-white noise
#'      assumption (instead of a weak-white noise assumption).
#'    \item **lowdiscrepancy**: Uses low-discrepancy sampling in the Monte-Carlo
#'      method. Note, low-discrepancy sampling will yield deterministic results.
#'  }
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#' hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#' 1-alpha quantile of the limiting distribution of the specified test's test statistic.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#'  limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#'  containing a short description of the test, the p-value, and additional information about the test if
#'  suppress_print_output = FALSE.
#'
#' @noRd
#' @keywords internal
#'
#' @references Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the autocovariance of a functional time series
#'  under conditional heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @examples
#' # b <- generate_brownian_motion(150)
#' # .multi_lag_test(b, lag = 5)
#' # .multi_lag_test(b, lag = 10, M = 50)
.multi_lag_test <- function(data, lag = 20, M=NULL,
                           method = c('iid', 'lowdiscrepancy'), alpha=0.05) {
  # Checks
  poss_methods <- c('iid', 'lowdiscrepancy')
  method <- .verify_input(method, poss_methods)

  # Setup
  data <- funts(data)
  f_data <- data$data

  N <- NCOL(f_data)
  J <- NROW(f_data)

  # Statistics
  if (method=='lowdiscrepancy') {
    # results <- V_WS_quantile(data$data, lag, alpha=alpha, M=M, low_disc=TRUE)

    bandwidth <- ceiling(0.25 * (N ^ (1/3)))

    covs <- autocovariance(center(funts(f_data))$data^2,
                                   lags = 1:lag, center=FALSE)
    mean_welch <- do.call("sum", covs) / J^2

    sum1 <- sum(
      sapply(1:lag, function(i,f_data,M,low_disc){
        MCint_eta_approx_i_j(f_data, i, i, M=M, low_disc=low_disc)
      }, f_data=f_data, M=M, low_disc=TRUE)
    )
    if (lag > 1) {
      for (i in 1:(lag-1)) {
        for (j in (i+1):lag) {
          if (abs(i-j) > bandwidth) next

          sum1 <- sum1 + 2 * MCint_eta_approx_i_j(f_data, i, j, M=M, low_disc=TRUE)
        }
      }
    }
    var_welch <- sum1


  } else if(method=='iid') {
    # results <- V_WS_quantile_iid(data$data, lag, alpha=alpha)

    cov <- autocovariance(f_data,0)

    mean_welch <- lag * (sum(diag(cov)) / J)^2
    var_welch <- lag * 2 * ( sum(cov^2) / J^2 )^2
  }

  ## Welch Parameters
  beta <- var_welch / (2 * mean_welch)
  nu <- 2 * mean_welch^2 / var_welch
  quantile <- beta * stats::qchisq(1 - alpha, nu)

  statistic <- sum(t_statistic_Q(f_data, 1:lag))
  p_val <- stats::pchisq(statistic / beta, nu, lower.tail = FALSE)

  # Return results
  list(statistic = statistic, quantile = quantile, p_value = p_val)
}


#' Spectral Density Test
#'
#' Computes the spectral hypothesis test under a user-specified kernel function and
#'  bandwidth; automatic bandwidth selection methods are provided.
#'
#' @description The "spectral" portmanteau test is based on the spectral density operator. It essentially measures
#'  the proximity of a functional time series to a white noise - the constant spectral density operator of an
#'  uncorrelated series. Unlike the "single-lag" and "multi-lag" tests, this test is not for general white noise
#'  series, and may not hold under functional conditionally heteroscedastic assumptions.
#'
#' @param data funts object or functional data matrix with observed functions in the columns
#' @param kernel A String specifying the kernel function to use. The currently supported kernels are the
#' 'Bartlett' and  'Parzen' kernels. The default kernel is 'Bartlett'.
#' @param bandwidth A String or positive Integer value which specifies the bandwidth to use. Currently admitted
#' string handles are 'static' which computes the bandwidth p via p = n^(1/(2q+1)) where n is the sample size
#' and q is the kernel order, or 'adaptive' which uses a bandwidth selection method that is based on the
#' functional data.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used for the test.
#'  The significance level is 0.05 by default. Note, the significance value is only ever used to compute the
#'  1-alpha quantile of the limiting distribution of the specified test's test statistic.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#'  limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#'  containing a short description of the test, the p-value, and additional information about the test if
#'  suppress_print_output = FALSE.
#'
#' @noRd
#' @keywords internal
#'
#' @references Characiejus V., & Rice G. (2019). A general white noise test based on kernel lag-window estimates of the
#' spectral density operator. Econometrics and Statistics, submitted.
#'
#' @references Chen W.W. & Deo R.S. (2004). Power transformations to induce normality and their applications.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66, 117-130.
#'
#' @examples
#' # b <- generate_brownian_motion(100)
#' # .spectral_test(b)
#' # .spectral_test(b, kernel = 'Parzen', bandwidth = 'adaptive')
#' # .spectral_test(b, kernel = 'Bartlett', bandwidth = 2)
.spectral_test <- function(data, kernel = 'Bartlett',
                          bandwidth = 'adaptive', alpha = 0.05) {
  ## TODO: Make kernels match rest of package
  data <- funts(data)
  f_data <- data$data

  J <- NROW(f_data)
  N <- NCOL(f_data)

  ## Get Kernel
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

  ## Get Bandwidth
  if (bandwidth == 'static') {
    bandwidth <- N^(1 / (2 * kernel_order + 1))
  } else if (bandwidth == 'adaptive') {
    bandwidth <- max(2, adaptive_bandwidth(f_data, kernel_string))
  } else if (!is.numeric(bandwidth)) {
    stop("Please see the documentation for valid bandwith arguments.")
  }

  ## Compute Statistic
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

  t_stat_term <- sigma_squared_hat^(-2) * C_hat_HS_norm[1] * sqrt(2 * D_n_k)
  untrans_num <- 2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq

  ### TODO: Add case for when H is not R
  beta <- 1 - (2/3) * sum(kernel_vals^2) * sum(kernel_vals^6) / (sum(kernel_vals^4)^2)
  t_stat <- ((2^(-1) * N * sigma_squared_hat^(-2) * spectral_distance_Q_sq)^beta -
               (C_n_k^beta + 2^(-1) * beta * (beta - 1) * C_n_k^(beta-2) * t_stat_term^2)) /
    (beta * C_n_k^(beta-1) * t_stat_term)

  # statistic <- spectral_t_statistic(data$data, kernel = kernel, bandwidth = bandwidth)

  list(statistic = t_stat,
       quantile = stats::qnorm(1 - alpha),
       p_value = 1 - stats::pnorm(t_stat),
       band = bandwidth)
}


#' Independence Test
#'
#' Computes the independence test with a user-specified number of principal components
#'  and range of lags.
#'
#' @details The "independence" portmanteau test is a test of independence and identical distribution based on a
#'  dimensionality reduction by projecting the data onto the most important functional principal components.
#'  It is based on the resulting lagged cross-variances. This test is not for general white noise series, and
#'  may not hold under functional conditionally heteroscedastic assumptions. Please consult the vignette for a
#'  deeper exposition, and consult the reference for a complete treatment.
#'
#' @param data funts object or functional data matrix with observed functions in the columns
#' @param components A positive Integer specifying the number of principal components to project the data on;
#'  ranked in order of importance (importance is determined by the proportion of the variance that is explained
#'  by the individual principal component.)
#' @param lag A positive Integer value, specifying the maximum lag to include - this can be seen as the bandwidth
#'  or lag-window.
#' @param alpha Numeric value between 0 and 1 specifying the significance level to be used in the specified
#'  hypothesis test. The default value is 0.05. Note, the significance value is only ever used to compute the
#'  1-alpha quantile of the limiting distribution of the specified test's test statistic.
#'
#' @return If suppress_raw_output = FALSE, a list containing the test statistic, the 1-alpha quantile of the
#'  limiting distribution, and the p-value computed from the specified hypothesis test. Also prints output
#'  containing a short description of the test, the p-value, and additional information about the test if
#'  suppress_print_output = FALSE.
#'
#' @noRd
#' @keywords internal
#'
#' @references Gabrys R., & Kokoszka P. (2007). Portmanteau Test of Independence for Functional Observations.
#' Journal of the American Statistical Association, 102:480, 1338-1348, DOI: 10.1198/016214507000001111.
#'
#' @examples
#' # b <- generate_brownian_motion(250)
#' # .independence_test(b, components = 3, lag = 5)
.independence_test <- function(data, components, lag, alpha = 0.05) {
  data <- funts(data)

  if ( (components < 1) | (components %% 1 != 0) ) {
    stop("The 'components parameter must be a positive integer.")
  }
  if ( (lag < 1) | (lag %% 1 != 0) ) {
    stop("The 'components lag must be a positive integer.")
  }

  N <- NCOL(data$data)
  J <- NROW(data$data)
  data <- center(data)

  suppressWarnings(
    pc_decomp <- ftsa::ftsm(rainbow::fts(1:J, data$data),
                            order = components, mean = FALSE)
    )
  scores <- pc_decomp$coeff
  C_0 <- crossprod(scores) / N
  c_h <- array(0, dim=c(components,components,lag))
  for (h in 1:lag) {
    for (k in 1:components) {
      for (l in 1:components) {
        score_uni <- 0
        for (t in 1:(N-h)) {
          score_uni <- score_uni + (scores[t,k] * scores[t+h,l])
        }
        c_h[k,l,h] <- score_uni / N
      }
    }
  }
  r_f_h <- r_b_h <- array(0, dim=c(components,components,lag))
  summand <- vector('numeric', lag)
  for (h in 1:lag) {
    r_f_h[,,h] <- solve(C_0) %*% c_h[,,h]
    r_b_h[,,h] <- c_h[,,h] %*% solve(C_0)
    summand[h] <- sum(r_f_h[,,h] * r_b_h[,,h])
  }
  Q_n <- N * sum(summand)
  p_val <- as.numeric(1 - stats::pchisq(Q_n, df = components^2 * lag))
  quantile <- as.numeric(stats::qchisq(1 - alpha, df = components^2 * lag))

  list(statistic = Q_n, quantile = quantile, p_value = p_val)
}


#' Imhof Test
#'
#' Returns the the SVD of the tensor c^hat_i_j(t,s,u,v) and the p-value computing
#'   the probability that the observed value of the statistic Q_h is larger than the 1-alpha
#'   quantile of the quadratic form in normal variables described in (15)
#'
#' @param data funts object or data that can be easily converted
#' @param lag the lag for which to compute the imhof test
#'
#' @return A list containing the SVD of tensor c^hat_i_j(t,s,u,v) and the p-value computing the
#'   probability that the observed value of the statistic Q_h is larger than the 1-alpha quantile
#'   of the quadratic form in normal variables.
#'
#' @noRd
#' @keywords internal
#'
#' @examples
#' # .imhof_test(funts(electricity[,1:50]),1)
.imhof_test <- function(data, lag) {
  data <- funts(data)

  if (!requireNamespace('tensorA')) {
    stop("Please install the 'tensorA' package to perform the imhof test.")
  }
  if (!requireNamespace('CompQuadForm')) {
    stop("Please install the 'CompQuadForm' package to perform the imhof test.")
  }
  if ( (lag < 1) | (lag %% 1 != 0) ) {
    stop("The 'lag' parameter must a positive integer.")
  }
  N <- NCOL(data$data)
  J <- NROW(data$data)

  t_statistic_val <- t_statistic_Q(data$data, lag)
  data <- center(data)

  tensor <- array(0, c(J, J, J, J))
  sum1 <- 0
  for (k in 1:(N-lag)) {
    tensor <- tensor + data$data[,k] %o% data$data[,k+lag] %o%
      data$data[,k] %o% data$data[,k+lag]
  }
  tensor <- tensor / N

  temp_tensor <- as.numeric(tensor)
  tensor_numeric <- tensorA::to.tensor(temp_tensor, c(J,J,J,J))
  names(tensor_numeric) <- c("a", "b", "c", "d")

  SVD <- tensorA::svd.tensor(tensor_numeric, i=c("a", "b"))
  eigenvalues <- as.numeric(SVD$d / (J^2))
  pval_imhof <- CompQuadForm::imhof(t_statistic_val, lambda = eigenvalues)$Qq

  list(statistic = t_statistic_val, p_value = pval_imhof)
}
