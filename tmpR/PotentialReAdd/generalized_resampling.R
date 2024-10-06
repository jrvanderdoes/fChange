#' Generalized Bootstrap and Permutation Methods Including Blocking
#'
#' Use this method for generalized resampling of a test statistic in
#'     determination of the proper cutoff value
#'
#' TODO:: Consider compiler::cmpfun(generalized_resampling)
#'
#' @param X A data frame of the evaluated fd data
#'     (rows are evaled points and columns are fd)
#' @param blockSize Numeric value indicating blocking sizes
#'
#' Value of 1 should be used for iid data while larger values can be used to
#'     account for dependence in the data
#'
#' @param fn Function that returns the test statistic value
#' @param iters Numeric value indicating number iterations used in bootstrapping
#' @param replace Boolean value indicating bootstrapping (T) or permutation (F)
#' @param alpha Numeric value in (0, 1) to get the correct percentile
#' @param silent Boolean to silence output helpful to a user
#' @param ... Optional arguments passed into \code{fn}
#'
#' @return Numeric value indicating the cutoff value
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup Data
#' data_KL <- generate_kl(
#'   ns = c(20, 20),
#'   eigsList = list(
#'     c(3, 2, 1, 0.5),
#'     c(3, 3, 2)
#'   ),
#'   basesList = list(
#'     fda::create.bspline.basis(nbasis = 4, norder = 4),
#'     fda::create.fourier.basis(nbasis = 2)
#'   ),
#'   meansList = c(0, 0.5),
#'   distsArray = c("Normal", "Binomial"),
#'   evals = seq(0, 1, 0.05),
#'   kappasArray = c(0, 0.5)
#' )
#' # Metric
#' #compute_Tn(data_KL$data,M =100) # Note value
#' #compute_Tn(data_KL$data,M =100) # Note different value
#' # Permutation Method for Tn (this will get 1-alpha quantile of iters)
#' generalized_resampling(
#'   X = data_KL$data,
#'   blockSize = 5,
#'   fn = compute_Tn, iters = 10, replace = F
#' )
#' }
generalized_resampling <- function(X, fn, blockSize=1, iters = 1000,
                                   replace = FALSE, alpha = 0.05,
                                   silent = FALSE, ...) {
  ## Get Function Value and estimate time
  st <- Sys.time()
  full_val <- fn(X, ...)
  en <- Sys.time()

  if (!silent) {
    cat(paste0(
      "Estimated time: ",
      round(difftime(en, st, "units" = "mins")[[1]] * iters, 2),
      " mins\n"
    ))
  }

  # ## Create Permuted Samples
  # n <- ncol(X)
  # idxGroups <- .getChunks(1:n, n / blockSize)
  # idxs <- sapply(1:iters, function(i, m, indxs, replace) {
  #   samps <- sample(1:m, replace = replace)
  #   unlist(indxs[samps], use.names = FALSE)
  # }, m = length(idxGroups), indxs = idxGroups, replace = replace)
  # # If it is a matrix/dataframe this already even rows
  # #     (will be with permutation test)
  # if (!(is.matrix(idxs) | is.data.frame(idxs))) {
  #   idxs <- .convertSamplesToDF(idxs)
  # }
  #
  # ## Sample via bootstrap
  # bssamples <- sapply(as.data.frame(idxs),
  #   function(loop_iter, fn, X1, ...) {
  #     fn(X1[, stats::na.omit(loop_iter)], ...)
  #   },
  #   X1 = X, fn = fn, ...
  # )

  ## Create Permuted Samples
  boot_dat <- .bootstrap(X,blockSize,iters,type,replace)
  bssamples <- unlist(lapply(boot_dat,fn, ...))


  list(
    "value" = full_val,
    "cutoff" = stats::quantile(bssamples, probs = c(1 - alpha))[[1]],
    "pval" = 1 - stats::ecdf(bssamples)(full_val),
    "BSSamples" = as.numeric(bssamples)
  )
}

