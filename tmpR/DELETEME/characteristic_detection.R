#' Characteristic Functions Change Point Detection - Simulation
#'
#' @param X Data.frame of the observations
#' @param M Numeric for the number of vectors to explore space
#' @param J Numeric for the reslolution to check
#' @param nSims Numeric for the number of simulations for distributions
#' @param h Numeric for bandwidth
#' @param K Function for bandwidth
#' @param space String for the space to explore space
#' @param silent Boolean to display outputs when running
#'
#' @return A list with the following elements:
#'    'Tn': List with ('pval') p-value, ('gamProcess') gamma process using in
#'        processes, and ('value') Tn statistic for data.
#'    'Mn': List with ('pval') p-value, ('gamProcess') gamma process using in
#'        processes, and ('value') Mn statistic for data.
#' @export
#'
#' @examples
#' characteristic_change_sim(electricity[,1:30],M = 10, J=25)
characteristic_change_sim <- function(X, M = 20, J=50,
                                      nSims = 1000, h = 3,
                                      K = bartlett_kernel, space = "BM",
                                      silent = FALSE) {
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_Mn <- compute_Mn(X = X, W = W, J = J)
  val_Tn <- compute_Tn(X = X, W = W, J = J)

  MJ <- M * J

  sqrtMat <- .compute_sqrtMat(X,W,J,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)
    # 0 0 0 0 (M times) ... .... 1 1 1 1 (M times)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    list(apply(results, MARGIN = 1, dot_integrate),
         apply(results, MARGIN = 1, max) )
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcessTn <- as.vector(unlist(gamProcess[1,]))
  gamProcessMn <- as.vector(unlist(gamProcess[2,]))

  list(
    'Tn'=list(
      "pval" = 1 - stats::ecdf(gamProcessTn)(val_Tn),
      "gamProcess" = gamProcessTn,
      "value" = val_Tn
    ),
    'Mn'=list(
      "pval" = 1 - stats::ecdf(gamProcessMn)(val_Mn$value),
      "gamProcess" = gamProcessMn,
      "value" = val_Mn$value
    ),
    'Location'=val_Mn$location
  )
}


#' Characteristic Functions Change Point Detection - Tn Simulation
#'
#' @inheritParams characteristic_change_sim
#'
#' @return List with three items:
#'  'pval': p-value of change point
#'  'gamProcess': Simulated gamma processes
#'  'value': Value of Tn test statistic for the data
#' @export
#'
#' @examples
#' characteristic_change_sim_Tn(electricity[,1:50],M=5,J=10,h=0)
characteristic_change_sim_Tn <- function(X,
                                     M = 20, J=50,
                                     nSims = 1000,
                                     h = 3,
                                     K = bartlett_kernel,
                                     space = "BM",
                                     silent = FALSE) {
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_Tn <- compute_Tn(X=X, W=W, J=J)
  #val_Tn <- compute_Tn_final1(X, W)
  #val_Tn <- compute_Tn_final2(X, W)

  MJ <- M * J

  sqrtMat <- .compute_sqrtMat(X,W,J,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    # todo / ncol(mvnorms1)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M,drop=FALSE])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    #apply(abs(gamVals)^2, MARGIN = 1, dot_integrate)
    apply(results, MARGIN = 1, dot_integrate)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Tn),
    "gamProcess" = gamProcess,
    "value" = val_Tn
  )
}


#' Characteristic Functions Change Point Detection - Mn Simulation
#'
#' @inheritParams characteristic_change_sim
#'
#' @return List with three items:
#'  'pval': p-value of change point
#'  'gamProcess': Simulated gamma processes
#'  'value': Value of Mn test statistic for the data
#' @export
#'
#' @examples
#' characteristic_change_sim_Mn(electricity[,1:50],M=5,J=10,h=0)
characteristic_change_sim_Mn <- function(X,
                                        M = 20, J=50,
                                        nSims = 1000,
                                        h = 3,
                                        K = bartlett_kernel,
                                        space = "BM",
                                        silent = FALSE) {
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_Mn <- compute_Mn(X = X, W = W, J=J)
  #val_Mn <- compute_Mn_final1(X, W)
  #val_Mn <- compute_Mn_final2(X, W)

  MJ <- M * J

  sqrtMat <- .compute_sqrtMat(X,W,J,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)
    # 0 0 0 0 (M times) ... .... 1 1 1 1 (M times)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    apply(results, MARGIN = 1, max)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Mn$value),
    "gamProcess" = gamProcess,
    "value" = val_Mn$value
  )
}


#############################################
#' Compute Characteristic Integrated Test Statistic
#'
#' Computes the integrated test statistic for the characteristic functional-based
#'  change point detection.
#'
#' @inheritParams characteristic_change_sim
#' @param W Vectors to explore the space against
#'
#' @return Numeric integrated test statistic
#' @export
#'
#' @examples
#' compute_Tn(electricity[,1:100])
compute_Tn <- function(X,
                       W=computeSpaceMeasuringVectors(M = 20, X = X, space = 'BM'),
                       J=50) {
  X <- funts(X)

  n <- ncol(X$data)

  Zn <- .Zn_final(W, X$data)
  ns <- .select_n(1:n, J)
  Zn <- Zn[ns,,drop=FALSE]
  # Integrate out W
  intVal <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  # Integrate observations
  dot_integrate(intVal)
}


#' Compute CUSUM value for integrated test statistic
#'
#' @inheritParams characteristic_change_sim
#'
#' @return Numeric CUSUM value
#'
#' @keywords internal
#' @noRd
.Zn_final <- function(W, X) {
  n <- ncol(X)
  fhat_vals <- as.matrix(.fhat_all(X, W))

  sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])
}



#' Compute Characteristic Maximized Test Statistic
#'
#' Computes the maximized test statistic for the characteristic functional-based
#'  change point detection.
#'
#' @inheritParams characteristic_change_sim
#' @param W Vectors to explore the space against
#'
#' @return List with the following
#'  \itemize{
#'    \item value: numeric maximized test statistic
#'    \item location: placement of maximized value
#'    \item allValues: test statistic value at each point
#'  }
#' @export
#'
#' @examples
#' compute_Mn(electricity[,1:100])
compute_Mn <- function(X,
                       W=computeSpaceMeasuringVectors(M = 20, X = X, space = 'BM'),
                       J=min(50,ncol(X))) {
  n <- ncol(X)

  Zn <- .Zn_final(W,X)
  ns <- .select_n(1:n,J)
  Zn <- Zn[ns,,drop=FALSE]
  # Integrate out W
  return_value <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  list(
    "value" = max(return_value),
    "location" = ns[which.max(return_value)],
    "allValues" = return_value
  )
}



#############################################

#' Compute the Square Root Matrix
#'
#' Computes the square root matrix for the characteristic functional-based
#'  change point detection
#'
#' @inheritParams characteristic_change_sim
#'
#' @return Matrix representing the square root matrix
#'
#' @keywords internal
#' @noRd
.compute_sqrtMat <- function(X,W,J,h,K){
  x <- seq(0,1,length.out=J)

  # Compute Gamma Matrix
  covMat <- .estimCovMat(X, W, x, h, K)

  # TODO:: Test "GauPro::sqrt_matrix()" which is slightly faster
  tryCatch({
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  }, error = function(e){
    covMat <- round(covMat, 15)
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  })

  Rfast::mat.mult(
    Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
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
#' @keywords internal
#' @noRd
.estimCovMat <- function(X, W, x,
                               h = 3, K = bartlett_kernel) {
  M <- ncol(W)

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
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
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
          lfun = funs[2][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
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
#' @keywords internal
#' @noRd
.estimD <- function(K, h, X, lfun, v, lfunp, vp) {
  iters <- (1 - ncol(X)):(ncol(X) - 1)

  # Move so pass in function values to avoid re-computation later
  fVals <- as.numeric(.estimf(X, lfun, v))
  fpVals <- as.numeric(.estimf(X, lfunp, vp))

  Kvals <- K(iters/h)
  data_tmp <- data.frame('K'=Kvals[Kvals>0],
                         'k'=iters[which(Kvals>0)])
  values <- apply(data_tmp, MARGIN = 1,
                  function(kInfo, X1, fVals, fpVals, mean1, mean2) {
    kInfo[1] * .estimGamma(
      k = kInfo[2], X = X1,
      fVals = fVals, fpVals = fpVals,
      mean1 = mean1, mean2 = mean2
    )
  },
  X1 = X, fVals = fVals, fpVals = fpVals,
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
#' @keywords internal
#' @noRd
.estimGamma <- function(k, X, fVals, fpVals, mean1, mean2) {

  if (k >= 0) {
    rs <- 1:(ncol(X) - k)
  } else {
    rs <- (1 - k):ncol(X)
  }

  sum(.estimR(rs, fVals, mean1) * Conj(.estimR(rs + k, fpVals, mean2)))
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
#' @keywords internal
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
#' @keywords internal
#' @noRd
.estimf <- function(Xr, lfun, v) {
  lfun((t(Xr) %*% v) / nrow(Xr))
}
