#' Compute Tn and Mn
#'
#' Function to compute both Tn and Mn.
#'
#' TODO:: Make this part of the original function
#'
#' @param X Data to compute the test statistics for
#' @param M Numeric for the number of vectors to measure space
#' @param W Data.frame of vectors to measure spaces
#' @param space String to indicate space for measure space
#' @param ... Additional details which are unused
#'
#' @return List with 'tn' for test statistic, 'mn' for test statistic,
#'  'location' for placement of change, and 'allmn' mn values at each
#'  point
#' @export
#'
#' @examples
#' compute_TnMn(electricity,M=5)
compute_TnMn <- function(X, M = 100000, W = NULL, space = "BM", ...) {
  # # Small bit of cpp Code for a function
  # Rcpp::cppFunction('ComplexMatrix col_cumsum(ComplexMatrix m) {
  #   for (int j = 0; j < m.ncol(); ++j) {
  #       for (int i = 1; i < m.nrow(); ++i) {
  #           m(i, j) = m(i, j) + m(i - 1, j);
  #       }
  #   }
  #   return m;
  # }')


  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }


  Zn2 <- abs(.Zn(W,X, n))^2
  intVal_tn <- 1 / M * sum(dot_integrate_col(Zn2))
  intVal_mn <- unname(unlist( 1 / M * rowSums(Zn2) ))

  list(
    "tn" = intVal_tn,
    "mn" = max(intVal_mn),
    "location" = which.max(intVal_mn),
    "allmn" = intVal_mn
  )
}


#' Resampling for Tn and Mn
#'
#' Compute resampling statistics for Tn and Mn
#'
#' TODO:: Make this part of the original function
#'
#' @param X data.frame of the data
#' @param blockSize Numeric for block length of the resampling
#' @param iters Numeric for the number of iterations for resampling
#' @param replace Boolean for replacement changing bootstrap to permutation
#' @param alpha Numeric between [0,1] for the significance level of alpha
#' @param silent Boolean for if anything should be output
#' @param ... Additional parameters
#'
#' @return List with the following items:
#'    'tn': Tn test statistic
#'    'mn': Mn test statistic
#'    'location': Estimate for the change location
#'    'cutoff_tn': Cutoff for the Tn test statistic
#'    'cutoff_mn': Cutoff for the Mn test statistic
#'    'pval_tn': p-value for Tn
#'    'pval_mn': p-value for Mn
#'    'data_bs': Bootstrap samples of the data
#' @export
#'
#' @examples
#' data_KL <- generate_kl(
#' ns = c(50, 50),
#' eigsList = list(
#'   c(3, 2, 1, 0.5),
#'   c(3, 3, 2)
#' ),
#' basesList = list(
#'   fda::create.bspline.basis(nbasis = 4, norder = 4),
#'   fda::create.fourier.basis(nbasis = 2)
#' ),
#' meansList = c(0, 0.5),
#' distsArray = c("Normal", "Binomial"),
#' evals = seq(0, 1, 0.05),
#' kappasArray = c(0, 0.5)
#' )
#' result <- TnMn_resampling(data_KL$data, M=1000)
TnMn_resampling <- function(X, blockSize = 1, iters = 1000,
                            replace = FALSE, alpha = 0.05, silent = FALSE, ...) {
  st <- Sys.time()
  full_val <- compute_TnMn(X, ...)
  en <- Sys.time()

  if (!silent) {
    cat(paste0(
      "Estimated time: ",
      round(difftime(en, st, "units" = "mins")[[1]] * iters, 2),
      " mins\n"
    ))
  }

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


  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter, X1, ...) {
                        compute_TnMn(X1[, stats::na.omit(loop_iter)], ...)
                      },
                      X1 = X, ...
                )

  co_tn <- stats::quantile(unname(unlist(bssamples[1,])), probs = c(1 - alpha))[[1]]
  co_mn <- stats::quantile(unname(unlist(bssamples[2,])), probs = c(1 - alpha))[[1]]

  list(
    "tn" = full_val$tn,
    "mn" = full_val$mn,
    "location" = full_val$location,
    "cutoff_tn" = co_tn,
    "cutoff_mn" = co_mn,
    "pval_tn" = 1 - stats::ecdf(unname(unlist(bssamples[1,])))(full_val$tn),
    "pval_mn" = 1 - stats::ecdf(unname(unlist(bssamples[2,])))(full_val$mn),
    "data_bs" = bssamples
  )
}
