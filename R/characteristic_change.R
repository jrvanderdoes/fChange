#' Characteristic Change Point Detection
#'
#' @inheritParams change
#'
#' @returns List with data like pvalue, test statistic,  location, simulation,
#'  and so forth.
#'
#' @noRd
#' @keywords internal
.change_characteristic <- function(X, statistic, critical,
                                   M = 20, J=50,
                                   nSims = 1000, h = 3,
                                   W = space_measuring_functions(X = X, M = 20, space='BM'),
                                   K = bartlett_kernel, #space = "BM",
                                   blocksize=1,
                                   resample_blocks = 'separate', replace = TRUE,
                                   alpha=0.05, ...) {
  X <- dfts(X)

  # Generate Noise
  if(is.null(W)){
    W <- space_measuring_functions(X = X, M = 20, ...)
  }
  M <- ncol(W)

  # Test Statistic
  tmp <- .characteristic_statistic(X = X$data,v=X$fparam,
                                    statistic=statistic, W = W, J = J,
                                    location = TRUE)
  stat <- tmp[1]
  location <- tmp[2]

  if(critical=='simulation'){

    sqrtMat <- .compute_sqrtMat(X$data,W,J,h,K)

    gamProcess <- c()
    nIters <- nSims / 100
    gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
      # (After trans + mult) Rows are iid MNV
      mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
      mvnorms <- Rfast::mat.mult(t(mvnorms), sqrtMat)

      gamVals <- mvnorms[,1:MJ] + complex(imaginary = 1) * mvnorms[,MJ + 1:MJ]

      # Integrate out M
      results <- matrix(NA,nrow=100, ncol=J)
      for(j in 1:J){
        results[,j] <- t(dot_integrate_col(t(abs(gamVals[,(j-1)*M + 1:M])^2)))
      }

      # Estimate value
      list(apply(results, MARGIN = 1, dot_integrate),
           apply(results, MARGIN = 1, max) )
    }, MJ = M * J, sqrtMat = sqrtMat, M=M, J=J)

    if(statistic=='Tn'){
      simulations <- as.vector(unlist(gamProcess[1,]))
    }else if(statistic=='Mn'){
      simulations <- as.vector(unlist(gamProcess[2,]))
    }

  } else if(critical=='resample'){

    simulations <- .bootstrap(X = X$data, blocksize = blocksize, M = nSims,
                         type = resample_blocks, replace = replace,
                         fn = .characteristic_statistic,
                         statistic=statistic, v=X$fparam, W = W, J = J)

  } else if(critical=='welch'){

    if(statistic != 'Tn') stop('Test statistic must be for Tn in Welch',call. = FALSE)

    ## Mean
    m_gamma <- 0
    for(i in 1:M){
      Cre <- .estimD(K=K, h=h, X=X$data,
                     lfun=cos, v=W[,i], lfunp=cos, vp=W[,i])
      Cim <- .estimD(K=K, h=h, X=X$data,
                     lfun=sin, v=W[,i], lfunp=sin, vp=W[,i])

      m_gamma <- m_gamma + ( Cre + Cim )
    }
    m_gamma <- m_gamma / (6*M)

    ## Variance
    sigma2_gamma <- sigma2_gamma1 <- 0
    for (i in 1:M) {
      val <- val1 <- 0
      for (j in 1:M){
        if (i==j) next
        Cre <- .estimD(K=K, h=h, X=X$data,
                       lfun=cos, v=W[,i], lfunp=cos, vp=W[,j])
        Cri <- .estimD(K=K, h=h, X=X$data,
                       lfun=cos, v=W[,i], lfunp=sin, vp=W[,j])
        Cim <- .estimD(K=K, h=h, X=X$data,
                       lfun=sin, v=W[,i], lfunp=sin, vp=W[,j])

        val <- val + (Cre^2+2*Cri^2+Cim^2)
      }
      val <- val / (M-1)
      sigma2_gamma <- sigma2_gamma + val
    }
    sigma2_gamma <- (2 / 90) * sigma2_gamma / M

    ## Welch Paramters
    beta <- sigma2_gamma / (2*m_gamma)
    nu <- (2*m_gamma^2) / sigma2_gamma

    val_cutoff <- beta * stats::qchisq(1 - alpha, df = nu)
    simulations <- beta * stats::rchisq(n = nSims, df = nu)

  }

  ## Setup and Return Data
  return_value <- list(
    "pvalue" = sum(stat <= simulations) / nSims,
    'location' = location#,
    # "statistic" = stat,
    # "simulations" = simulations
  )
  # if(critical=='welch') {
  #   return_value[['critical']] <- val_cutoff
  #   return_value[['alpha']] <- alpha
  # }

  return_value
}


#' Compute Characteristic Test Statistic
#'
#' Computes the integrated/maximum test statistic for characteristic
#'  functional-based change point detection.
#'
#' @inheritParams characteristic_change_sim
#' @param W Vectors to explore the space against
#' @param v Intraday evaluation points
#' @param location Boolean if the location should also be return. Default is FALSE.
#'
#' @return Numeric test statistic and potentially the location
#'
#' @noRd
#' @keywords internal
.characteristic_statistic <- function(
    X, statistic='Tn', v=seq(0,1,length.out=nrow(X)),
    W = space_measuring_functions(M = 20, X = X, space = 'BM'), J = 50,
    location = FALSE, all.stats=FALSE){

  n <- ncol(X)

  # Compute Zn
  fhat_vals <- as.matrix(.fhat_all(X, W))
  Zn <- sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])

  ns <- .select_n(1:n, J)
  Zn <- Zn[ns,,drop=FALSE]

  # Integrate out W
  noW <- dot_integrate_col(t(abs(Zn)^2))

  if(statistic=='Tn'){
    statistic <- dot_integrate(noW)

  } else if (statistic=='Mn'){
    statistic <- max(noW)
  }

  if(!location) return(statistic)

  c( statistic, ns[which.max(noW)] )
}


#' Get fhat Estimates
#'
#' This computes the values in fHat, i.e. cutoff for each possible value.
#'  Note, a starting 0 is omitted. See usage in .Zn.
#'
#' @inheritParams .characteristic_statistic
#'
#' @return Vector of numerics for \eqn{\hat{f}(v,x)}
#'
#' @keywords internal
#' @noRd
.fhat_all <- function(X, W) {
  dot_col_cumsum(exp(1 / nrow(X) * complex(imaginary = 1) * t(X) %*% as.matrix(W)) / ncol(X))
}


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
