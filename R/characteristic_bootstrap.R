#' Title
#'
#' @param X
#' @param statistic
#' @param M
#' @param J
#' @param space
#' @param blockSize
#' @param iters
#' @param replace
#' @param alpha
#' @param silent
#'
#' @return
#' @export
#'
#' @examples
detect_changepoint_bootstrap <- function(X, statistic=c('Tn','Mn'),
                                         M=20, J=50, space='BM',
                                         blockSize=1, iters = 1000,
                                         replace = FALSE, alpha = 0.05,
                                         silent = FALSE) {
  # Test Statistics
  if(length(statistic)!=1){
    stop('Choose "Tn" or "Mn" as the test statistic')
  }else if(statistic=='Tn'){
    fn <- compute_Tn_final
  }else if(statistic=='Mn'){
    fn <- compute_Mn_final
  }else{
    stop('Choose "Tn" or "Mn" as the test statistic')
  }

  # Noise
  W <- computeSpaceMeasuringVectors(M,space,X)

  ## Get Function Value and estimate time
  st <- Sys.time()
  full_val <- fn(X, W=W, J=J)[[1]]
  en <- Sys.time()

  if (!silent) {
    cat(paste0(
      "Estimated time: ",
      round(difftime(en, st, "units" = "mins")[[1]] * iters, 2),
      " mins\n"
    ))
  }

  ## Create Permuted Samples
  n <- ncol(X)
  idxGroups <- .getChunks(1:n, n / blockSize)
  idxs <- sapply(1:iters, function(i, m, indxs, replace) {
    samps <- sample(1:m, replace = replace)
    unlist(indxs[samps], use.names = FALSE)
  }, m = length(idxGroups), indxs = idxGroups, replace = replace)
  # If it is a matrix/dataframe this already even rows
  #     (will be with permutation test)
  if (!(is.matrix(idxs) | is.data.frame(idxs))) {
    idxs <- .convertSamplesToDF(idxs)
  }

  ## Sample via bootstrap
  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter, fn, X1, W, J) {
                        fn(X1[, stats::na.omit(loop_iter)], W=W, J=J)[[1]]
                      },
                      X1 = X, fn = fn, W=W, J=J
  )

  list(
    "value" = full_val,
    "cutoff" = stats::quantile(bssamples, probs = c(1 - alpha))[[1]],
    "pval" = 1 - stats::ecdf(bssamples)(full_val),
    "BSSamples" = as.numeric(bssamples)
  )
}
