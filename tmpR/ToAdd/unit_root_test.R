#' Functional Unit Root Test
#'
#' Unit root test for functional time series with different methods on determining the critical values of the test statistic. The Monte Carlo method was constructed in Chen and Pun (2021), while the bootstrap-based methods have not been validated in the literature (although such an option is provided, please use them at your own risk).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @param ALPHA Significance level of the test. The default value is 5\%.
#' @param TYPE Type of hypothesis test. The default value 'TYPE="R"' represents the right-tailed test, which is used when the alternative hypothesis is trend stationarity model. 'TYPE="L"' represents the left-tailed test, which is used when the alternative hypothesis is simple stationarity model or AR(1) model.
#' @param METHOD Method to determine the critical value of the test statistic. The default value 'METHOD="MC"' represents the Monte Carlo method. 'METHOD="SB"' represents the simple bootstrap method and 'METHOD="MBB"' represents the moving block bootstrap method.
#' @param K Kernel function in the estimation of the long-run covariance function, which is only effective in the Monte Carlo method. The default function is 'default_kernel' function in this package.
#' @param h_power Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth, which is only effective in the Monte Carlo method. The default value is 2/5.
#' @param est_ev Number of the largest eigenvalues chosen to estimate the limiting distribution, which is only effective in the Monte Carlo method. The default value is the sample size 'N'.
#' @param MCNsim Number of Monte Carlo datasets generated in the Monte Carlo method, which is only effective in Monte Carlo method. The default value is 10000.
#' @param bm_set A vector of independent simulated data generated from the function 'dataset_bm', which is only effective and essential in Monte Carlo method.
#' @param M Number of bootstrap datasets generated in the bootstrap method, which is only effective in bootstrap methods. The default value is 1000.
#' @param b Block length used in the moving block bootstrap method, which is only effective in the moving block bootstrap method. The default value is ceiling((2N)^(1/3)), where 'N' is the sample size.
#'
#' @return The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of functional random walk.
#' @export
#'
#' @references Chen, Y., & Pun, C. S. (2021). Functional Unit Root Test.
#'  Available at SSRN.
#'
#' @examples
#' res <- compute_unit_root_test(
#'   generate_brownian_motion(100, v=seq(0,1,length.out=20)))
#' res <- compute_unit_root_test(
#'   generate_brownian_motion(100, v=seq(0,1,length.out=20)), statistic='Mn')
#' res1 <- compute_unit_root_test(
#'   generate_brownian_motion(1000, v=seq(0,1,length.out=30)))
#' res2 <- compute_unit_root_test(electricity)
compute_unit_root_test <- function(X, type='trend', method="MC", boot_method='seperate',
                                   M=2500, h = 3, K=bartlett_kernel,
                                   TVE=0.95, replace=TRUE){
  stop('Under Work')
  # https://github.com/cran/STFTS/blob/master/R/furt.R
  # https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3761262#:~:text=In%20this%20paper,%20we%20propose%20a%20unit%20root%20test%20for
  # method="MC"; boot_method='seperate'
  # M=1000; h = 3; TVE=0.95; replace=TRUE
  # test='trend'#c('trend','unitroot','stationarity')
  # X <- generate_brownian_motion(5000,v=seq(0,1,length.out=20))
  # X <- electricity

  X <- .check_data(X)
  N <- ncol(X$data)
  r <- nrow(X$data)

  hat_eta <- diff(X)$data

  .compute_unit_root_statistic  <- function(hat_eta,N,r){

    # 1/N^2 * sum(dot_integrate_col( (X$data-X$data[,1])^2 ) )

    ZN <- 1/sqrt(N) * cumsum(funts(hat_eta))$data
    # TODO: uneven cols integration
    dot_integrate(dot_integrate_col(ZN^2))
  }

  stat <- .compute_unit_root_statistic(hat_eta,N,r)

  if (method=="MC"){
    pca_X <- pca(funts(hat_eta), TVE = TVE)

    sim_stats <- sapply(1:M,function(m, eigs, v){
      sum(eigs *
            dot_integrate_col(generate_brownian_motion(length(eigs), v = v)$data^2) )
    },eigs=pca_X$sdev^2, v=seq(0,1,0.01))
    # cat(paste0(round(TN,4),' - ',
    #            round(head(sim_stats)[1],4),', ', round(head(sim_stats)[2],4),
    #            ' - ', mean(TN<sim_stats)))
  } else if(method=="BS"){
    boot_X <- .bootstrap(hat_eta, blockSize=h, M=M,
                         type=boot_method, replace=replace)
    sim_stats <- rep(NA, M)
    for(i in 1:M){
      sim_stats[i] <- .compute_unit_root_statistic(funts(boot_X[[i]]),N,r)
    }
  }

  if(type=='trend'){
    list('statistic'=stat, 'pvalue'=mean(stat < sim_stats) )
  }else if (type=='unitroot' || type=='stationarity'){
    list('statistic'=stat, 'pvalue'=mean(stat > sim_stats) )
  }
}
