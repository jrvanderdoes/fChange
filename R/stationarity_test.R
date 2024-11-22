#' Functional Stationarity Test
#'
#' Stationarity test for functional time series with different methods on
#'  determining the critical values of the test statistic. The Monte Carlo
#'  method was constructed in Horvath et al. (2014), while the bootstrap-based
#'  methods have not been validated in the literature (although such an option
#'  is provided, please use them at your own risk).
#'
#' @param X The functional time series being tested, inputted in a matrix form
#'  with each row representing each observation of the functional data values
#'  on equidistant points of any pre-specified interval.
#' @param statistic description
#' @param method description
#' @param boot_method description
#' @param M description
#' @param h description
#' @param TVE description
#' @param replace description
#'
#' @return The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of local stationarity.
#' @export
#'
#' @references Horvath, L., Kokoszka, P., & Rice, G. (2014). Testing
#'  stationarity of functional time series. Journal of Econometrics,
#'  179(1), 66-82.
#'
#' @examples
#' res <- compute_stationarity_test(
#'   generate_brownian_motion(100,v=seq(0,1,length.out=20)))
#' res <- compute_stationarity_test(
#'   generate_brownian_motion(100,v=seq(0,1,length.out=20)),statistic='Mn')
#' res1 <- compute_stationarity_test(
#'   generate_brownian_motion(1000,v=seq(0,1,length.out=30)))
#' res2 <- compute_stationarity_test(electricity)
compute_stationarity_test <-
  function(X, statistic = 'Tn', method="MC", boot_method='seperate',
           M=1000, h = 3, TVE=0.95, replace=TRUE){
  X <- .check_data(X)

  N <- ncol(X$data)
  r <- nrow(X$data)

  stat <- .compute_stationary_test_stat(X,statistic)

  if(method=="MC"){
    pca_X <- pca(X, TVE = TVE)

    if(statistic == 'Tn'){
      # TODO: Uneven integration for Tn
      sim_stats <- sapply(1:M,function(m, eigs, v){
        sum(eigs *
              dot_integrate_col(generate_brownian_bridge(length(eigs), v = v)$data^2) )
      },eigs=pca_X$sdev^2, v=X$intraobs)
    }else if(statistic=='Mn'){
      sim_stats <- sapply(1:M,function(m, eigs, v){
        bbs <- generate_brownian_bridge(length(eigs), v = v)$data
        sum(eigs * dot_integrate_col(
          (bbs - matrix(dot_integrate_col(bbs), nrow = length(v),
                        ncol = length(eigs), byrow = TRUE))^2 ) )
      },eigs=pca_X$sdev^2, v=X$intraobs)
    }else{
      stop('Statistic is incorrect',call. = FALSE)
    }

  } else if(method=="BS"){
    boot_X <- .bootstrap(X, blockSize=h, M=M,
                         type=boot_method, replace=replace)
    sim_stats <- rep(NA, M)
    for(i in 1:M){
      sim_stats[i] <- .compute_stationary_test_stat(funts(boot_X[[i]]),statistic)
    }
  }

  list('statistic' = stat, 'pvalue' = sum(stat <= sim_stats)/M)
}


#' Compute Stationarity Test Statistic
#'
#' Internal function to compute stationarity test statistic
#'
#' @inheritParams compute_stationarity_test
#'
#' @return Numeric test statistic
#'
#' @keywords internal
#' @noRd
.compute_stationary_test_stat <- function(X,statistic){
  X <- .check_data(X)
  N <- ncol(X$data)
  r <- nrow(X$data)

  # Test Statistics
  Zn <- 1/sqrt(N) * ( cumsum(X)$data -
                        matrix((1:N)/N, nrow = r, ncol = N, byrow = TRUE) * rowSums(X$data) )
  if(statistic == 'Tn'){
    # TODO: Uneven integration for Tn
    stat <- dot_integrate(dot_integrate_col(Zn^2))
  }else if(statistic=='Mn'){
    # TODO:: Fix uneven integration
    stat <- dot_integrate_uneven(dot_integrate_col(
      t(Zn - matrix(dot_integrate_col(t(Zn)),nrow=r,ncol=N,byrow = FALSE) )^2),
      r = X$intraobs )
  }else{
    stop('Statistic is incorrect',call. = FALSE)
  }

  stat
}
