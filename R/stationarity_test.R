#' Functional Stationarity Test
#'
#' Stationarity test for functional time series with different methods on
#'  determining the critical values of the test statistic. The Monte Carlo
#'  method was constructed in Horvath et al. (2014), while the resample-based
#'  methods have not been validated in the literature (use the provided option
#'  at your discretion).
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param statistic String for test statistic. Options are integrated (\code{Tn})
#'  and supremum (\code{Mn}). The default is \code{Tn}.
#' @param critical String for method of determining the critical values. Options
#'  are \code{simulation} and \code{resample}. Default is \code{simulation}.
#' @param perm_method String for method of resampling. Options are \code{separate}
#'  for block resampling and \code{overlapping} for sliding window.
#'  Default is \code{separate}.
#' @param M Numeric for number of simulation to use in determining the null
#'  distribution. Default is 1000.
#' @param blocksize Numeric for blocksize in resample test. Default is \eqn{2N^{1/5}}.
#' @param TVE Numeric for total variance explained when using PCA for
#'  eigenvalues. Default is 1.
#' @param replace Boolean if replacement should be used for resample test. Thus,
#'  this defines if a bootstrapped or permuted test is used. Default is TRUE.
#'
#' @return List with the following elements:
#'  \enumerate{
#'    \item pvalue: p-value for the stationarity test.
#'    \item statistic: test statistic from the test.
#'    \item simulations: simulations which define the null distribution.
#'  }
#' @export
#'
#' @references Horvath, L., Kokoszka, P., & Rice, G. (2014). Testing
#'  stationarity of functional time series. Journal of Econometrics,
#'  179(1), 66-82.
#'
#' @examples
#' res <- stationarity_test(
#'   generate_brownian_motion(100,v=seq(0,1,length.out=20)),
#'   critical='resample', statistic='Mn')
#' res2 <- stationarity_test(electricity)
stationarity_test <-
  function(X, statistic = 'Tn', critical=c('simulation','resample'),
           perm_method='separate',  M=1000, blocksize = 2*ncol(X)^(1/5),
           TVE=1, replace=TRUE){
    critical <- .verify_input(critical, c('simulation','resample'))

    X <- dfts(X)

  N <- ncol(X$data)
  r <- nrow(X$data)

  stat <- .compute_stationary_test_stat(X$data, X$fparam, statistic)

  if(critical=="simulation"){
    pca_X <- pca(X, TVE = TVE)

    if(statistic == 'Tn'){
      sim_stats <- sapply(1:M,function(m, eigs, v){
        sum(eigs *
              dot_integrate_col(
                generate_brownian_bridge(length(eigs), v = v)$data^2, v) )
      },eigs=pca_X$sdev^2, v=X$fparam)
    }else if(statistic=='Mn'){
      sim_stats <- sapply(1:M,function(m, eigs, v){
        vv <- v
        bbs <- generate_brownian_bridge(length(eigs), v = vv)$data
        sum(eigs * dot_integrate_col(
          (bbs - matrix(dot_integrate_col(bbs), nrow = length(vv),
                        ncol = length(eigs), byrow = TRUE))^2, vv ) )
      },eigs=pca_X$sdev^2, v=X$fparam)
  }else{
      stop('Statistic is incorrect',call. = FALSE)
    }

  } else if(critical=="resample"){
    sim_stats <- .bootstrap(X$data, blocksize=blocksize, M=M,
                            type=perm_method, replace=replace,
                            fn=.compute_stationary_test_stat,
                            statistic=statistic, v=X$fparam)
  }

  list('pvalue' = sum(stat <= sim_stats)/M,
       'statistic' = stat,
       'simulations' = sim_stats)
}


#' Compute Stationarity Test Statistic
#'
#' Internal function to compute stationarity test statistic
#'
#' @inheritParams stationarity_test
#'
#' @return Numeric test statistic
#'
#' @keywords internal
#' @noRd
.compute_stationary_test_stat <- function(X, v, statistic){
  N <- ncol(X)
  r <- nrow(X)

  # Test Statistics
  ## TODO:: CUsum with data.frame
  Zn <- 1/sqrt(N) * ( cumsum(dfts(X))$data -
                        matrix((1:N)/N, nrow = r, ncol = N, byrow = TRUE) * rowSums(X) )
  if(statistic == 'Tn'){
    stat <- dot_integrate(dot_integrate_col(Zn^2, v))
  }else if(statistic=='Mn'){
    stat <- dot_integrate(
      dot_integrate_col( ( Zn -
                           matrix(dot_integrate_col(t(Zn)), nrow = length(v),
                                  ncol = N, byrow = FALSE))^2, v) )
  }else{
    stop('Statistic is incorrect',call. = FALSE)
  }

  stat
}
