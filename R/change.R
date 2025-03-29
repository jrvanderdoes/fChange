#' Change Point Detection
#'
#' Change point detection for dfts objects. Various change point methods
#'  are given, where single or multiple changes can be detected. Multiple change
#'  extensions currently include binary segmentation and elbow plots.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param method Method to compute change point. Options include:
#'  'characteristic', 'mean', 'robustmean', 'eigenjoint', 'eigensingle', 'trace',
#'  'covariance', 'projmean', and 'projdistribution'.
#' @param statistic String for the test statistic type: integrated, \code{Tn},
#'  or supremum, \code{Mn}.
#' @param critical String for method of computing threshold. Options are
#'    'simulation', 'resample', and 'welch'. Not all ways to compute the critical
#'    thresholds are implemented for every method.
#' @param type String for the type of change point detection, single change
#'  ('single'), binary segmentation ('segmentation'), or elbow plots ('elbow').
#' @param resample_blocks String indicating the type of resample test to use.
#'  Using \code{separate} gives blocks which are separate while \code{overlapping}
#'  creates overlapping or sliding windows. When \code{blocksize=1} then these
#'  will be identical.
#' @param replace Boolean for using a permutation or bootstrapped statistic when
#'  \code{critical='resample'}.
#' @param max_changes Integer as the max number of changes to search when using
#'  type is \code{elbow}.
#' @param changes Vector of change points to be given to the eigen test if the data
#'  should be centered on these values first.
#' @param blocksize Integer for the width of the blocks when using a resampling
#'  test. Can use [adaptive_bandwidth()] if additional guidance is desired.
#' @param eigen_number Which eigenvalue or the number of eigenvalues which should be checked
#'  in the eigenvalue tests.
#' @param h Number of lags used when computing long run covariance estimates. Used in
#'  mean, characteristic, and eigenvalue tests.
#' @param M Number of simulations or permutations for critical values
#' @param J Resolution (J) in the characteristic method. The number of vectors
#'  is defined by \code{W}.
#' @param W Space measuring functions used in characteristic method to explore
#'  the functional space.
#' @param K Kernel function for use in characteristic, mean, eigen, covariance
#'  and projmean.
#' @param alpha Significance in \[0,1\] for Welch approximation.
#' @param cov.res Resolution to use when computing covariance kernel changes.
#' @param weighting Weights used in covariance kernel method and pcadistribution.
#' @param TVE Total variance explained for projmean and projdistribution.
#' @param trim_function Trimming to be used in multiple change methods.
#' @param errors Type of errors used in elbow plot. Options are L2 and Trace.
#' @param recommendation_change_points Number of lags forward to examine in deciding automated
#'  elbow plot recommendation.
#' @param recommendation_improvement Significant drop to look for in deciding automated elbow
#'  plot recommendation.
#' @param silent.binary Boolean if output should be printed when running binary
#'  segmentation.
#'
#' @returns
#' When type is single, returns a list:
#'  \enumerate{
#'    \item pvalue: p-value for detection of a change point.
#'    \item location: location of the most likely change.
#'  }
#'  When type is elbow:
#'  \enumerate{
#'    \item information: data.frame with the information on each change and the
#'      decrease in variability.
#'    \item plots: list of plots showing the variability decrease or improvement
#'    \item suggestion: list with plot and algorithmic change suggestion. The
#'      suggested changes are also returned.
#'  }
#'  When type is segmentation a data.frame with the locations and p-values is
#'  returned.
#'
#' @export
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2018). Detecting and dating structural
#'  breaks in functional data without dimension reduction. Journal of the Royal
#'  Statistical Society. Series B, Statistical Methodology, 80(3), 509-529.
#'
#' @references Wegner, L., Wendler, M. Robust change-point detection for
#'  functional time series based on U-statistics and dependent wild bootstrap.
#'  Stat Papers (2024).
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2020). Structural break
#'  analysis for spectrum and trace of covariance operators. Environmetrics
#'  (London, Ont.), 31(1)
#'
#' @references Horvath, L., Rice, G., & Zhao, Y. (2022). Change point analysis
#'  of covariance functions: A weighted cumulative sum approach. Journal of
#'  Multivariate Analysis, 189, 104877-.
#'
#' @references Berkes, I., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Detecting changes in the mean of functional observations}
#'  Journal of the Royal Statistical Society, Series B 71, 927-946
#'
#' @references Aue, A., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Estimation of a change-point in the mean function of functional data}
#'  Journal of Multivariate Analysis 100, 2254-2269.
#'
#' @references Huskova, M., & Meintanis, S.G. (2006). Change Point Analysis
#'  based on Empirical Characteristic Functions. Metrika, 63, 145-168.
#'
#' @examples
#' res <- fchange(electricity$data[,1:20],method='characteristic',critical = 'welch')
fchange <- function(X,
                   method=c('characteristic','mean','robustmean','eigenjoint',
                            'eigensingle','trace',
                            'covariance','projmean','projdistribution'),
                   statistic=c('Tn','Mn'),
                   critical=c('simulation','resample','welch'),
                   type=c('single','segmentation','elbow'),
                   resample_blocks = 'separate', replace=TRUE,
                   max_changes=min(ncol(X),20),
                   changes=NULL,
                   blocksize = 2*ncol(X)^(1/5),
                   eigen_number = 3, h = 2*ncol(X)^(1/5),
                   M = 1000, J=50,
                   W = space_measuring_functions(X = X, M = 20, space='BM'),
                   K = bartlett_kernel,
                   alpha=0.05, cov.res = 30, weighting = 1/4, TVE=0.95,
                   trim_function = function(X) { 0 },
                   errors='L2', recommendation_change_points = 2,
                   recommendation_improvement = 0.15,
                   silent.binary = FALSE){

  # Check Data
  X <- dfts(X)
  fake_return <- list('pvalue'=1,'location'=NA)
  if(dim(X)[2]==1) return(fake_return)
  method <- .verify_input(method,
                          c('characteristic','mean','robustmean','eigenjoint',
                            'eigensingle','trace',
                            'covariance','projmean','projdistribution'))
  statistic <- .verify_input(statistic, c('Tn','Mn'))
  critical <- .verify_input(critical, c('simulation','resample','welch'))
  type <- .verify_input(type, c('single','segmentation','elbow'))
  max_changes <- round(max_changes)

  # Detect Changes
  if(type=='single'){
    result <- switch(method,
                     characteristic={
                       # TODO:: Add if multiple test statistics given..
                       .change_characteristic(
                         X = X, statistic=statistic, critical=critical,
                         J=J,
                         nSims = M, h = h,
                         K = K, W = W,
                         blocksize=blocksize, resample_blocks = resample_blocks,
                         replace = replace, alpha=alpha)
                     },
                     mean={
                       # TODO:: Welch Approximation
                       .change_mean(data = X, statistic=statistic,
                                              critical=critical, M = M, h = h,
                                              K = K, blocksize=blocksize,
                                              type = resample_blocks, replace = replace)
                     },
                     robustmean={
                       # TODO:: Welch check
                       # TODO:: bandwidth
                       .change_robust(X, statistic = statistic,
                                                bandwidth = NA,
                                                m = M, threshold = critical)
                     },
                     eigenjoint={
                       # TODO:: Welch check
                       .change_eigen(X = X, d = eigen_number, h = h,
                                     changes = changes,
                                     statistic = statistic,
                                     test='joint',
                                     critical = critical,
                                     blocksize = blocksize,
                                     M = M,K = K,
                                     type = resample_blocks,
                                     replace = replace)
                     },
                     eigensingle={
                       # TODO:: Welch check
                       .change_eigen(X = X, d = eigen_number, h = h,
                                     changes = changes,
                                     statistic = statistic,
                                     test='individual',
                                     critical = critical,
                                     blocksize = blocksize,
                                     M = M,K = K,
                                     type = resample_blocks,
                                     replace = replace)
                     },
                     trace={
                       # TODO:: Welch check
                       .change_trace(X = X, changes = changes, M = M,
                                     statistic = statistic,
                                     critical = critical,
                                     blocksize = blocksize,
                                     replace = replace,
                                     type = resample_blocks)
                     },
                     covariance={
                       # TODO:: Welch check
                       # TODO:: Mn
                       .change_covariance_kernel(X=X, statistic=statistic,
                                                 critical=critical,
                                                 kappa = weighting, len = cov.res,
                                                 blocksize=blocksize, M=M,
                                                 resample_blocks=resample_blocks,
                                                 replace=replace,
                                                 K=K)
                     },
                     projmean={
                       # TODO:: Welch check
                       .change_pca_mean(X=X, statistic=statistic,
                                        critical=critical,
                                        TVE=TVE, M=M, K=K,
                                        blocksize=blocksize,
                                        resample_blocks=resample_blocks,
                                        replace=replace )
                     },
                     projdistribution={
                       if(critical != 'resample')
                         stop('Only resample setup for this method currently',call. = FALSE)
                       .change_pca_distribution(X=X, statistic=statistic, critical=critical,
                                                TVE = TVE, gam = weighting, M = M,
                                                blocksize = blocksize, resample_blocks = resample_blocks,
                                                replace = replace)
                     },
                     {
                       # Default
                     }
    )
  } else if(type=='segmentation'){
    result <-
      .binary_segmentation(X=X, method=method,
                           statistic=statistic, critical=critical,
                           resample_blocks = resample_blocks, replace=replace,
                           changes=changes, blocksize = blocksize,
                           eigen_number=eigen_number, h=h, M = M, J=J, W = W, K = K,
                           alpha=alpha, cov.res = cov.res, weighting = weighting,
                           TVE=TVE, trim_function = trim_function,
                           silent = silent.binary)


  } else if(type=='elbow'){

    result <-
      .elbow_method(X, method=method, W=W,
                   trim_function = trim_function,
                   max_changes = max_changes,
                   errors = errors,
                   K=K, d=eigen_number, h=h, weighting=weighting,
                   recommendation_change_points = recommendation_change_points,
                   recommendation_improvement = recommendation_improvement,
                   TVE=TVE)
  }

  result
}
