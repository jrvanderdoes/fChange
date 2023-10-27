#' Detect Change Point
#'
#' This method detects changes in the data through simulating the null
#'  distribution. To do so, the covariance matrix is estimated based on the
#'  data, against some noise vectors. After estimation, a sample from the null
#'  distribution can be simulated. The creation of the covariance matrix and
#'  simulating of the null sample is completed many times, estimating the
#'  null distribution. The test statistic is compared to this. This method is
#'  computationally intensive and we found no reason to prefer it over
#'  `detect_changepoint_singleCov()`.
#'
#' @param X Numeric data.frame of functional data observations--rows for
#'    evaluated values and columns indicating FD
#' @param nSims (Optional) Integer indicating the number of realizations of
#'     Gaussian processes to compute. Default is 100
#' @param x (Optional) Vector of locations of observations. Default is equally
#'     spaced observations on (0, 1) with the same number of observations as FDs.
#' @param h (Optional) Integer indicating amount of lag to consider. Default is
#'     3.
#' @param K (Optional) Function for the kernel function to use. Default is the
#'     barlett_kernel.
#' @param space (Optional) String indicating the space used in the spanning vectors.
#'     Default is 'BM'.
#' @param TN_M (Optional) Integer indicating the number of iteractions used to compute
#'     T_N for the real data. Default is 10000.
#' @param Cov_M (Optional) Integer indicating the number of vectors used to create
#'     each value in the covariance matrix. Default is 25.
#' @param silent (Optional) Boolean indicating it the output should be supressed.
#'     Default is FALSE.
#'
#' @return List with three entries:
#'  \enumerate{
#'    \item pval: pvalue based on the data and estimated Gaussian processes
#'    \item gamProcess: Vector of estimated test statistics based on data
#'    \item value: Test statistic for data
#'  }
#' @export
#'
#' @examples
#' cp_res <- detect_changepoint(
#'   electricity[, 1:8],
#'   nSims = 1, x = seq(0, 1, length.out = 5),
#'   h = 0, K = bartlett_kernel
#' )
detect_changepoint <- function(X, nSims = 500, x = seq(0, 1, length.out = 40),
                               h = 3, K = bartlett_kernel, space = "BM",
                               TN_M = 10000, Cov_M = 75, silent = FALSE) {
  # Setup Vars
  gamProcess <- rep(NA, nSims)
  value <- compute_Tn(X, M = TN_M, space = space)

  MJ <- Cov_M * length(x)

  for (i in 1:nSims) {
    if (!silent) cat(paste0(i, ", "))
    # Compute Gamma Matrix
    W <- computeSpaceMeasuringVectors(Cov_M, space, X)
    covMat <- .estimCovMat(X = X, x = x, M = Cov_M, h = h, K = K, W = W)

    mvnorms <- suppressWarnings(mvtnorm::rmvnorm(1, sigma = covMat))

    gamVals <- mvnorms[1, 1:MJ] + complex(imaginary = 1) * mvnorms[1, MJ + 1:MJ]

    # Estimate value
    gamProcess[i] <- .approx_int(abs(gamVals)^2)
  }

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(value),
    "gamProcess" = gamProcess,
    "value" = value
  )
}

#' Estimate null and detect change point based on a single covariance matrix
#'
#' This method detects changes in the data through simulating the null
#'  distribution. To do so, the covariance matrix is estimated based on the
#'  data, against some noise vectors. After estimation, samples from the null
#'  distribution can be simulated. These samples can be combine to estimate the
#'  null distribution. The test statistic is then compared to this. This method
#'  has been shown to be effective and is far faster than
#'  `detect_changepoint()`.
#'
#' @inheritParams detect_changepoint
#'
#' @return List with three entries:
#'  \enumerate{
#'    \item pval: pvalue based on the data and estimated Gaussian processes
#'    \item gamProcess: Vector of estimated test statistics based on data
#'    \item value: Test statistic for data
#'  }
#' @export
#'
#' @examples
#' cp_res <- detect_changepoint_singleCov(
#'   electricity[, 1:10],
#'   nSims = 100, x = seq(0, 1, length.out = 3),
#'   h = 0, K = bartlett_kernel, silent = FALSE
#' )
detect_changepoint_singleCov <- function(X, nSims = 2000, x = seq(0, 1, length.out = 40),
                                         h = 3, K = bartlett_kernel, space = "BM",
                                         silent = FALSE, TN_M = 100000, Cov_M = 75) {
  # Source Cpp File to speed up matrix computations
  # Rcpp::sourceCpp("R/matrixMult.cpp")

  # Determine Number of Iterations
  val_Tn <- compute_Tn(X, M = TN_M, space = space)

  # Generate Noise
  W <- computeSpaceMeasuringVectors(Cov_M, space, X)

  # Variables
  MJ <- Cov_M * length(x)

  # Compute Gamma Matrix
  covMat <- .estimCovMat(X, x, Cov_M, h, K, W)
  covMat <- round(covMat, 15)
  covMat_svd <- La.svd(covMat)
  sqrtD <- diag(sqrt(covMat_svd$d))
  sqrtD[is.na(sqrtD)] <- 0
  # sqrtMat <- eigenMapMatMult(
  #   eigenMapMatMult(covMat_svd$u, sqrtD),t(covMat_svd$v))
  sqrtMat <- Rfast::mat.mult(
    Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
  )
  rm(W, covMat, covMat_svd, sqrtD)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, lx) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- sapply(1:100, function(m, x) {
      stats::rnorm(x)
    }, x = 2 * MJ)
    mvnorms <- Rfast::mat.mult(sqrtMat, mvnorms)

    gamVals <- mvnorms[1:MJ, ] + complex(imaginary = 1) * mvnorms[MJ + 1:MJ, ]
    # gamVals <- mvnorms[(M+1):MJ,] + complex(imaginary = 1)*mvnorms[MJ+1:(MJ-M),]

    # Estimate value
    apply(abs(t(gamVals))^2, MARGIN = 1, .approx_int)
  }, MJ = MJ, sqrtMat = sqrtMat, lx = length(x))
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Tn),
    "gamProcess" = gamProcess,
    "value" = val_Tn
  )
}


#' Estimate the Covariance Matrix
#'
#' This (internal) function computes the covariance matrix of some FD data.
#'
#' @inheritParams detect_changepoint
#' @param M (Optional) Integer indicating the number of vectors used to create
#'     each value in the covariance matrix. Default is 25.
#' @param W (Optional) Data.frame of vectors for measuring space. If NULL then
#'     then a Guassian measure is generated. Default is NULL.
#'
#' @return Data.frame for covariance based on Gaussian measure and given data X
#'
#' @noRd
.estimCovMat <- function(X, x = seq(0, 1, length.out = nrow(X)),
                         M = 25, h = 3, K = bartlett_kernel, W = NULL) {
  # Setup random vectors
  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M, "BM", X)
  }
  if (ncol(W) != M) stop("Number of Vectors not M")

  D11 <- D12_tD21 <- D22 <- data.frame(matrix(ncol = M, nrow = M))
  funs <- c(cos, sin)
  for (i in 1:M) {
    for (j in i:M) {
      D11[i, j] <- D11[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[1][[1]], vp = W[, j]
        )
      D12_tD21[i, j] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i], lfunp = funs[2][[1]], vp = W[, j]
        )
      if (j > i) {
        D12_tD21[j, i] <- .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, j],
          lfunp = funs[2][[1]], vp = W[, i]
        )
      }

      D22[i, j] <- D22[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[2][[1]], v = W[, i], lfunp = funs[2][[1]], vp = W[, j]
        )
    }
  }

  tmp1 <- x %*% t(x)
  minDF <- matrix(1, nrow = length(x), ncol = length(x))
  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      minDF[i, j] <- (min(x[i], x[j]) - tmp1[i, j])
    }
  }

  # Mults value from minDF to D11, then next value of minDF to D11, ...
  # rbind(cbind(fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D11)),
  #             fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D12))),
  #       cbind(fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D21)),
  #             fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D22))))
  fastmatrix::kronecker.prod(
    as.matrix(minDF),
    rbind(
      cbind(as.matrix(D11), as.matrix(D12_tD21)),
      cbind(as.matrix(t(D12_tD21)), as.matrix(D22))
    )
  )
}


#' Estimate Long-run Covariance (D) matrix
#'
#' This (internal) function
#'
#' @inheritParams detect_changepoint
#' @param K Function for the kernel function to use.
#' @param h Integer indicating amount of lag to consider.
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#' @param lfunp Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param vp Numeric vector from L2 space with same number of observations as
#'     FD observation.
#'
#' @return Data.frame with autocovariance value based on data given
#'
#' @noRd
.estimD <- function(K, h, X, lfun, v, lfunp, vp) {
  iters <- (1 - ncol(X)):(ncol(X) - 1)

  # Move so pass in function values to avoid re-computation later
  fVals <- as.numeric(.estimf(X, lfun, v))
  fpVals <- as.numeric(.estimf(X, lfunp, vp))

  values <- sapply(iters, function(k, K, h, X1, fVals, fpVals, mean1, mean2) {
    Kval <- K(k, h)

    ifelse(Kval == 0,
      0,
      Kval * .estimGamma(
        k = k, X = X1,
        fVals = fVals, fpVals = fpVals,
        mean1 = mean1, mean2 = mean2
      )
    )
  },
  K = K, h = h, X1 = X, fVals = fVals, fpVals = fpVals,
  mean1 = mean(fVals), mean2 = mean(fpVals)
  )

  sum(values) / ncol(X)
}


#' Estimate Autocovariance (gamma) Function
#'
#' This (internal) function computes the autocovariance (gamma) function.
#'
#' @inheritParams .estimD
#' @param k Integer indicating FD object to consider
#'
#' @return Numeric autocovariance value for given data
#'
#' @noRd
.estimGamma <- function(k, X, fVals, fpVals, mean1, mean2) {
  # tmp <- 0

  if (k >= 0) {
    rs <- 1:(ncol(X) - k)
  } else {
    rs <- (1 - k):ncol(X)
  }

  sum(.estimR(rs, fVals, mean1) * .estimR(rs + k, fpVals, mean2))
}


#' Estimate R-hat
#'
#' This (internal) function estimates R-hat, that is lfun(<X_r,v>) - Y where Y
#'     demeans the data, typically mean(lfun(<X,v>)).
#'
#' @inheritParams .estimGamma
#' @param r Integer indicating the FD object of interest in the data X
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#' @param meanVal (Optional) Value indicating the meanVal. This can be is
#'     computed beforehand for speed. Default of NA computes mean(f-hat)
#'
#' @return Numeric value of R-hat for fd object of interest
#'
#' @noRd
.estimR <- function(r, fVals, meanVal = NA) {
  if (is.na(meanVal)) meanVal <- mean(fVals)

  fVals[r] - meanVal
}


#' Estimate f-hat
#'
#' This (internal) function computes estimate of f-hat, that is function(<X,v>).
#'
#' @param Xr Numeric vector/data.frame indicating FD observation(s).
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#'
#' @return Numeric value(s) indicating function value for each FD object given
#'
#' @noRd
.estimf <- function(Xr, lfun, v) {
  lfun((t(Xr) %*% v) / nrow(Xr))
}
