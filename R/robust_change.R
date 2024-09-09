#' Robust Change Point Detection
#'
#' @param n length of observation
#' @param d dimension of observation
#' @param ht heavy tailed scenario (1) or normal-distribution (0)
#' @param m number of bootstrap repetitions
#' @param S number of simulation runs
#' @param quantiles desired quantiles
#'
#' @return A list with the following elements
#' \itemize{
#'  \item Wilcox_pvalue: p-value for the Wilcox test statistic
#'  \item Cusum_pvalue: p-value for the cusum test statistic
#'  \item Wilcox_bandwidth: Bandwidth for Wilcox test
#'  \item Cusum_bandwidth: Bandwidth for cusum test
#' }
#' @export
#'
#' @references Wegner, L., Wendler, M. Robust change-point detection for
#'  functional time series based on U-statistics and dependent wild bootstrap.
#'  Stat Papers (2024). https://doi.org/10.1007/s00362-024-01577-7
#'
#' @examples
#' robust_change(n = 200, d = 100, ht = 0,  m = 100, S = 10,  quantiles = 0.95)
robust_change <- function(X, m){#, type=c('Wilcox','Cusum')){

  # Create Proper Data and Bandwidths
  Obs <- X$data
  # if(!is.na(pmatch('Wilcox',type))){
  #
  # }
  Obs_tilde_h <- make_Obs_tilde_h(Obs)
  Obs_tilde_C <- make_Obs_tilde_C(Obs)
  bw_h <- adaptive_bw(Obs_tilde_h)
  bw_C <- adaptive_bw(Obs_tilde_C)

  # calculate quadratic spectral cov. matrices
  ker_h <- kernel(bw_h,ncol(Obs))
  A_h <- toeplitz(ker_h)
  Re_A_h <- getRealSQM(A_h)
  Re_A_C <- Re_A_h

  # Calculate own Kernel for CUSUM, as needed
  if(bw_C != bw_h){
    ker_C <- kernel(bw_C,ncol(Obs))
    A_C <- toeplitz(ker_C)
    Re_A_C <- getRealSQM(A_C)
  }

  # find max_k U_{n,k}^(t) for t=1,...,m
  hC_Obs <- make_hC_Obs(Obs)
  ## TODO:: fill_U is slow
  U_hC <- fill_U(A_h=Re_A_h,
                 A_C=Re_A_C,
                 norm=rnorm(m*ncol(Obs)),
                 hC_Obs=hC_Obs,
                 m=m, n=ncol(Obs), d=nrow(Obs) )
  maxis <- find_max(U_hC, m, ncol(Obs))

  # find max_k U_{n,k}
  Te <- fill_T(hC_Obs, ncol(Obs), nrow(Obs))
  T_max <- apply(Te,MARGIN = 2, max)

  list("Wilcox_pvalue"=1-ecdf(maxis[,1])(T_max[1]),
       "Cusum_pvalue"=1-ecdf(maxis[,2])(T_max[2]),
       "Wilcox_bandwidth"=bw_h,
       "Cusum_bandwidth"=bw_C)
}


