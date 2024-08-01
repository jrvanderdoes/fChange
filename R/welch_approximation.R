#' Detect the Changepoint Based on Welch Approximation
#'
#' @param X Data.frame of the data
#' @param M Numeric for the number of vectors to measure space
#' @param J Numeric for the level of discretion
#' @param h Numeric for bandwidth
#' @param K Function for bandwidth
#' @param space String for the space to measure data
#' @param alpha Numeric for significance
#'
#' @return List with the following items:
#'    'cutoff': Cutoff calculated for Tn using Welch approximation,
#'    'value': Test statistic Tn for the data
#'    'detected': Boolean if Tn is significant or not
#'    'alpha': Significance of the result
#' @export
#'
#' @examples
#' detect_changepoint_Welch(electricity)
detect_changepoint_Welch <- function(X,
                                     M = 20, J=50,
                                     h = ncol(X)^(1 / 3),
                                     K = bartlett_kernel,
                                     space = "BM",
                                     alpha = 0.05){
  # Generate Noise
  W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  W1 <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  val_cutoff <- compute_Welch(X, alpha = alpha, W = W, W1 = W1,
                          M = M, h = h, K = K )
  val_Tn <- compute_Tn_final(X, W, J)

  list(
    "cutoff" = val_cutoff,
    "value" = val_Tn,
    'detected' = val_Tn > val_cutoff,
    'alpha'=alpha
  )
}

#' Welch Approximation to Tn
#'
#' This function approximates the Tn statistic using the Welch approximation.
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param alpha (Optional) Numeric value in (0, 1) for significance. Default is 0.05
#' @param TVal (Optional) Numeric value indicate the number of FDs. Default is the
#'     number of columns in X.
#' @param W (Optional) Data.frame for the functions to integrate against
#'     used for both muHat and sigmaHat2. Default of NULL will use Brownian
#'     motion.
#' @param W1 (Optional)  data.frame for the functions to integrate against
#'     used for sigmaHat2. Default of NULL will use Brownian motion.
#' @param M (Optional) Numeric value for number of functions to integrate against
#'     (ncol(W)). Default is 100.
#' @param h (Optional) Numeric value indicating the amount of blocking. Default
#'     is Tval^(1/3)
#' @param K (Optional) Function used for CHat, indicating Kernel. Default is the
#'     Bartlett Kernel.
#' @param ... Unused parameters, but added for use in other functions.
#'
#' @return Numeric value indicating cutoff from Welsh approximation
#' @export
#'
#' @examples
#' compute_Welch(electricity[,1:40], h = 0)
compute_Welch <- function(X, alpha = 0.05,
                          W = NULL, W1 = NULL, M = 100,
                          h = ncol(X)^(1 / 3),
                          K = bartlett_kernel ) {
  # Setup Vars
  TVal = ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M, "BM", X)
  }
  if (is.null(W1)) {
    W1 <- computeSpaceMeasuringVectors(M, "BM", X)
  }

  muVals <- sigma2Vals <- rep(0, M)

  for (i in 1:M) {
    v <- W[, i]
    v1 <- W1[, i]

    CHat <- .computeLRVEstimate(X = X, v1 = v, v2 = v, TVal = TVal, h = h, K = K)
    CHat1 <- .computeLRVEstimate(X = X, v1 = v, v2 = v1, TVal = TVal, h = h, K = K)

    muVals[i] <- CHat
    sigma2Vals[i] <- CHat1 * Conj(CHat1) # pairs
  }

  muHat <- 1 / 6 * mean(muVals)
  sigma2Hat <- (1 / 90) * 2 * mean(sigma2Vals)

  betaHat <- Re(sigma2Hat / (2 * muHat))
  nuHat <- Re((2 * muHat^2) / sigma2Hat)

  betaHat * stats::qchisq(1 - alpha, df = nuHat)
}


#' Compute Long Run Variance
#'
#' This (internal) function computes long run variance estimates
#'
#' @param h Numeric value indicating the amount of blocking.
#' @param TVal Numeric value indicate the number of FDs.
#' @param gammaHat Vector of numerics for gamma estimates at each lag.
#' @param K Function for the Kernel.
#'
#' @return Numeric estimating LRV value
#'
#' @noRd
.computeLRVEstimate <- function(X, v1, v2, TVal, h, K) {
  iters <- (1 - TVal):(TVal - 1)

  # Move so pass in function values to avoid re-computation later
  eps1 <- exp(1 / nrow(X) * complex(imaginary = 1) * (t(X) %*% v1))
  eps2 <- exp(1 / nrow(X) * complex(imaginary = 1) * (t(X) %*% v2))

  values <- sapply(iters, function(k, K, h, X1, TVal, eps1, eps2, mean1, mean2) {
    Kval <- K(k, h)

    ifelse(Kval == 0,
      0,
      Kval * .estimWelchGamma(
        k = k, X = X1, TVal = TVal,
        eps1 = eps1, eps2 = eps2,
        mean1 = mean1, mean2 = mean2
      )
    )
  },
  K = K, h = h, X1 = X, TVal = TVal, eps1 = eps1, eps2 = eps2,
  mean1 = mean(eps1), mean2 = mean(eps2)
  )

  sum(values)
}

#' Estimate Autocovariance Estimates
#'
#' This (internal) function estimates autocovariance at single point
#'
#' @param eps1 Complex vector indicating exp(i, <X,v>)
#' @param eps2 Complex vector indicating exp(i, <X,v2>)
#' @param l Integer indicating lags
#' @param TVal Numeric value indicate the number of FDs.
#'
#' @return Numeric indicating autocovariance value at a particular lag
#'
#' @noRd
.estimWelchGamma <- function(k, X, TVal, eps1, eps2, mean1, mean2) {
  if (k >= 0) {
    rs <- (k + 1):TVal
  } else {
    rs <- (1 - k):TVal
  }

  1 / (TVal) * sum((eps1[rs - abs(k)] - mean1) * Conj(eps2[rs] - mean2))
}
