#' Functional Hypothesis Tests for Functional Data
#'
#' Computes a variety of portmanteau hypothesis tests for functional data in the
#'  form of dfts objects.
#'
#' @details The "single"-lag portmanteau test assesses the significance of
#'  empirical lagged autocovariance operators at a single lag \code{lag}.
#'  It tests the null hypothesis that the lag-h autocovariance operator is equal
#'  to 0. The test is designed for stationary functional time-series, and is
#'  valid under conditional heteroscedasticity conditions.
#'
#' The "multi"-lag portmanteau test assesses the cumulative significance of
#'  empirical lagged autocovariance operators, up to a user-selected maximum lag
#'  \code{lag}. It tests the null hypothesis that the first lag-h autocovariance
#'  operators, \eqn{h=1,\dots,lag}, is equal to 0. The test is designed for
#'  stationary functional time-series, and is valid under conditional
#'  heteroscedasticity conditions.
#'
#' The "spectral" portmanteau test measures the proximity of a functional time
#'  series to a white noise. Comparison is made to the constant spectral
#'  density operator of an uncorrelated series. The test is not for general
#'  white noise series, and may not hold under functional conditionally
#'  heteroscedastic assumptions.
#'
#' The "independence" portmanteau test measures independence and identical
#'  distribution based lagged cross-variances from dimension reduction using
#'  functional principal components analysis. The test is not for general white
#'  noise series, and may not hold under functional conditionally
#'  heteroscedastic assumptions.
#'
#' The "imhof" portmanteau test is an analogue of the "single-lag" test. While
#'  the "single-lag" test computes the limiting distribution of the test
#'  statistic via a Welch-Satterthwaite approximation, the "imhof" test directly
#'  computes the coefficients of the quadratic form in normal variables. Hence,
#'  the test is computationally expensive.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param test A String specifying the hypothesis test. Currently available
#'  tests are: 'variety', 'single-lag', 'multi-lag', 'spectral', 'independence',
#'  and 'imhof'.
#' @param lag A positive integer to specify the lag, or maximum lag, of
#'  interest. Only used for the "single-lag", "multi-lag", "independence", and
#'  "imhof" tests.
#' @param M Numeric to specify the number of Monte-Carlo or resampled
#'  simulations to use for the limiting distributions.
#' @param method String indicating the method for the \code{single} test,
#'  options include:
#'  \itemize{
#'    \item **iid**: The
#'      hypothesis test will use a strong-white noise assumption (instead of a
#'      weak-white noise assumption).
#'    \item **resample**: The hypothesis
#'      test is evaluated by approximating the limiting distribution of the test
#'      statistic via a (block) resampling process.
#'  }
#'  Additional methods are forthcoming.
#' @param kernel Kernel function for spectral test or estimation of covariance.
#' @param block_size Numeric to specify block size for resampling tests.
#' @param bandwidth Numeric for bandwidth of covariance estimation. If left null,
#'  with be defined by \eqn{N^{1 / (2 * \text{KO} + 1)}} where KO is the order
#'  of the selected kernel.
#' @param components Number of functional principal components to use in the
#'  independence test.
#' @param resample_blocks String indicating the type of resample test to use.
#'  Using \code{separate} gives blocks which are separate while \code{overlapping}
#'  creates overlapping or sliding windows. When \code{blocksize=1} then these
#'  will be identical.
#' @param replace Boolean for using a permutation or bootstrapped statistic when
#'  \code{method='resample'}.
#' @param alpha Numeric value for significance in \[0,1\].
#'
#' @return List with results dependent on the test. In general, returns the pvalue,
#'  statistic, and simulations/quantile.
#' @export
#'
#' @references Kim, M., Kokoszka P., & Rice G. (2023) White noise testing for
#'  functional time series. Statist. Surv., 17, 119-168, DOI: 10.1214/23-SS143
#'
#' @references Characiejus, V., & Rice, G. (2020). A general white noise test
#'  based on kernel lag-window estimates of the spectral density operator.
#'  Econometrics and Statistics, 13, 175–196.
#'
#' @references Kokoszka P., & Rice G., & Shang H.L. (2017). Inference for the
#'  autocovariance of a functional time series under conditional
#'  heteroscedasticity. Journal of Multivariate Analysis, 162, 32-50.
#'
#' @references Zhang X. (2016). White noise testing and model diagnostic
#'  checking for functional time series. Journal of Econometrics, 194, 76-95.
#'
#' @references Gabrys R., & Kokoszka P. (2007). Portmanteau Test of Independence
#'  for Functional Observations. Journal of the American Statistical Association,
#'  102:480, 1338-1348, DOI: 10.1198/016214507000001111.
#'
#' @references Chen W.W. & Deo R.S. (2004). Power transformations to induce
#'  normality and their applications. Journal of the Royal Statistical Society:
#'  Series B (Statistical Methodology), 66, 117-130.
#'
#' @examples
#' b <- generate_brownian_motion(250)
#' res0 <- portmanteau_tests(b, test = 'single', lag = 10)
#' res1 <- portmanteau_tests(b, test = 'multi', lag = 10, alpha = 0.01)
#' res2 <- portmanteau_tests(b, test = 'spectral', kernel = bartlett_kernel,
#'                           bandwidth = NULL, alpha = 0.05)
#' res3 <- portmanteau_tests(b, test = 'spectral', alpha = 0.1,
#'                           kernel = parzen_kernel,
#'                           bandwidth = adaptive_bandwidth(b, kernel=parzen_kernel))
#' res4 <- portmanteau_tests(b, test = 'independence', components = 3, lag = 3)
#' res5 <- portmanteau_tests(b, test = 'single', lag = 1, M = 250)
portmanteau_tests <- function(X, test = c('variety', 'single', 'multi',
                                      'spectral', 'independence', 'imhof'),
                       lag=5, M=1000,
                       method = c('iid','bootstrap'),
                       kernel = bartlett_kernel,
                       block_size = NULL, bandwidth=NULL,
                       components = 3,
                       resample_blocks='separate', replace = FALSE, alpha=0.05) {
  X <- dfts(X)
  poss_tests <- c('variety', 'single', 'multi', 'spectral',
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

  ## Define bandwidth for passing
  if (!is.null(bandwidth) ) {
    bandwidth <- bandwidth + 0
  }else if(test %in% c('variety','single','spectral')){
    ## Get Kernel
    tmp <- .kernel_details(kernel,name=deparse(substitute(kernel)))
    kernel_order <- tmp$order

    bandwidth <- ncol(X)^(1 / (2 * kernel_order + 1))
  }

  ## RUN TESTS
  if (test == 'variety') {
    m <- as.table(matrix(0, ncol = 1, 10))
    colnames(m) <- c('pvalue')
    rownames(m) <- c('single-lag, lag = 1', 'single-lag, lag = 2',
                     'single-lag, lag = 3', 'multi-lag, lag = 5',
                     'multi-lag, lag = 10', 'multi-lag, lag = 20',
                     'spectral, static bandwidth',
                     'spectral, adaptive bandwidth',
                     'independence, 3 components, lag = 3',
                     'independence, 16 components, lag = 10')

    m[1] <- portmanteau_tests(X, test = 'single', lag = 1, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue
    m[2] <- portmanteau_tests(X, test = 'single', lag = 2, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue
    m[3] <- portmanteau_tests(X, test = 'single', lag = 3, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue

    m[4] <- portmanteau_tests(X, test = 'multi', lag = 5, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue
    m[5] <- portmanteau_tests(X, test = 'multi', lag = 10, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue
    m[6] <- portmanteau_tests(X, test = 'multi', lag = 20, alpha=alpha,
                              M=M, method=method, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue

    m[7] <- portmanteau_tests(X, test = 'spectral', lag=lag, method=method,
                              alpha=alpha, M=M, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue
    m[8] <- portmanteau_tests(X, test = 'spectral', lag=lag, method=method,
                              alpha=alpha, M=M, kernel=kernel,
                              block_size = block_size, bandwidth=bandwidth,
                              components = components,  replace =  replace,
                              resample_blocks=resample_blocks)$pvalue

    m[9] <- portmanteau_tests(X, test = 'independence', method=method,
                       components = 3, lag = 3, alpha=alpha,
                       M=M, kernel=kernel,
                       block_size = block_size, bandwidth=bandwidth,
                       replace =  replace,
                       resample_blocks=resample_blocks)$pvalue

    m[10] <- portmanteau_tests(X, test = 'independence', method=method,
                        components = 16, lag = 10, alpha=alpha,
                        M=M, kernel=kernel,
                        block_size = block_size, bandwidth=bandwidth,
                        replace =  replace,
                        resample_blocks=resample_blocks)$pvalue
    results <- m
  } else if (test == 'multi') {
    method <- 'iid'
    results <- .multi_lag_test(X, lag, M=M, method=method, alpha=alpha)
  } else if (test == 'single') {
    results <- .single_lag_test(X, lag, alpha=alpha,
                    M=M, method=method,
                    block_size=block_size,
                    resample_blocks = resample_blocks,
                    kernel=kernel, replace=replace)
  } else if (test == 'spectral') {
    results <- .spectral_test(X, kernel = kernel, bandwidth = bandwidth,
                              alpha = alpha)
  } else if (test == 'independence') {
    results <- .independence_test(X, components = components, lag = lag, alpha = alpha)
  } else if (test == 'imhof') {
    input <- readline("The imhof test is computationally expensive. \n
                      Press [enter] if you would like to continue.")
    if (input != '') {
      stop("User cancelled the test.")
    }
    results <- .imhof_test(X, lag)
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
#' @inheritParams portmanteau_tests
#' @param method string indicating type of test.
#'  \itemize{
#'    \item **iid**: The hypothesis test will use a strong-white
#'      noise assumption (instead of a weak-white noise assumption).
#'    \item **lowdiscrepancy**: The hypothesis test will usea
#'      low-discrepancy sampling in the Monte-Carlo method. Note,
#'      low-discrepancy sampling will yield deterministic results. (BROKEN)
#'    \item **bootstrap**: The hypothesis test is done by approximating the
#'      limiting distribution of the test statistic via a block bootstrap
#'      process.
#'  }
#' @param M Positive integer value. Number of Monte-Carlo simulations for the Welch-Satterthwaite approximation.
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
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .single_lag_test(generate_brownian_motion(100), lag = 1)
#'  }
.single_lag_test <- function(
    data, lag=1, alpha=0.05, method = c('iid', 'lowdiscrepancy', 'bootstrap'),
    M=500, resample_blocks = 'separate',
    kernel=bartlett_kernel, replace=FALSE,
    block_size=adaptive_bandwidth(data$data, kernel = kernel) ){

  data <- dfts(data)

  poss_methods <- c('iid', 'lowdiscrepancy', 'bootstrap')
  method <- .verify_input(method, poss_methods)

  if(method=='bootstrap') {
    stats_distr <- .bootstrap(data, blocksize = block_size,
                                    M = M, type =  resample_blocks,
                                    replace = replace, lag=lag, fn=t_statistic_Q)

    quantile <- stats::quantile(as.numeric(stats_distr), 1 - alpha)
    statistic <- t_statistic_Q(data$data, lag)
    pvalue <- sum(statistic <= stats_distr) / length(stats_distr)

    results <- list(pvalue = as.numeric(pvalue),
                    statistic = as.numeric(statistic),
                    quantile = as.numeric(quantile),
                    block_size = block_size)
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
#' @param data A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
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
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .multi_lag_test(generate_brownian_motion(100), lag = 5)
#'  }
.multi_lag_test <- function(data, lag = 20, M=NULL,
                           method = c('iid', 'lowdiscrepancy'), alpha=0.05) {

  # Checks
  poss_methods <- c('iid', 'lowdiscrepancy')
  method <- .verify_input(method, poss_methods)

  # Setup
  data <- dfts(data)
  f_data <- data$data

  N <- NCOL(f_data)
  J <- NROW(f_data)

  # Statistics
  if (method=='lowdiscrepancy') {
    # results <- V_WS_quantile(data$data, lag, alpha=alpha, M=M, low_disc=TRUE)

    bandwidth <- ceiling(0.25 * (N ^ (1/3)))

    covs <- autocovariance(center(dfts(f_data))$data^2,
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
  list(pvalue = p_val, statistic = statistic, quantile = quantile)
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
#' @param data A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
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
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .spectral_test(generate_brownian_motion(100))
#'  }
.spectral_test <- function(data, kernel = bartlett_kernel,
                          bandwidth = NULL, alpha = 0.05, order=NULL) {

  J <- NROW(data$data)
  N <- NCOL(data$data)

  ## Get Bandwidth
  if (is.null(bandwidth) ) {
    bandwidth <- N^(1 / (2 * order + 1))
  }

  ## Compute Statistic
  data_inner_prod <- crossprod(data$data) / J
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

  list(
    pvalue = 1 - stats::pnorm(t_stat),
    statistic = t_stat,
    quantile = stats::qnorm(1 - alpha),
    bandwidth = bandwidth)
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
#' @param data dfts object or functional data matrix with observed functions in the columns
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
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .independence_test(generate_brownian_motion(100), components=3, lag=5)
#'  }
.independence_test <- function(data, components, lag, alpha = 0.05) {
  data <- dfts(data)

  if ( (components < 1) | (components %% 1 != 0) ) {
    stop("The 'components parameter must be a positive integer.")
  }
  if ( (lag < 1) | (lag %% 1 != 0) ) {
    stop("The 'components lag must be a positive integer.")
  }

  N <- NCOL(data$data)
  J <- NROW(data$data)
  data <- center(data)

  # TODO:: Convert to internal PCA
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

  list(pvalue = p_val, statistic = Q_n, quantile = quantile)
}


#' Imhof Test
#'
#' Returns the the SVD of the tensor c^hat_i_j(t,s,u,v) and the p-value computing
#'   the probability that the observed value of the statistic Q_h is larger than the 1-alpha
#'   quantile of the quadratic form in normal variables described in (15)
#'
#' @param data A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param lag the lag for which to compute the imhof test
#'
#' @return A list containing the SVD of tensor c^hat_i_j(t,s,u,v) and the p-value computing the
#'   probability that the observed value of the statistic Q_h is larger than the 1-alpha quantile
#'   of the quadratic form in normal variables.
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .imhof_test(dfts(electricity$data[,1:50]),1)
#'  }
.imhof_test <- function(data, lag) {
  data <- dfts(data)

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
  eigenvalues <- as.numeric(SVD$d / J^2)
  pval_imhof <- CompQuadForm::imhof(t_statistic_val, lambda = eigenvalues)$Qq

  list(pvalue = pval_imhof, statistic = t_statistic_val)
}
