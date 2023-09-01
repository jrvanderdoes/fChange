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
#' # Note, TVal generally bigger than 1.
#' welch_approximation(electricity, TVal = 1, h = 0)
welch_approximation <- function(X, alpha = 0.05, TVal = ncol(X),
                                W = NULL, W1 = NULL, M = 100,
                                h = floor(TVal^(1 / 3)), K = bartlett_kernel, ...) {
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
  sigma2Hat <- (1 / 90) * 2 * mean(sigma2Vals) / length(sigma2Vals)

  betaHat <- Re(sigma2Hat / (2 * muHat))
  nuHat <- Re((2 * muHat^2) / sigma2Hat)

  betaHat * stats::qchisq(1 - alpha, df = nuHat) # / nrow(X)
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

  1/rs * sum( (eps1[rs - abs(k)] - mean1) * Conj(eps2[rs]-mean2) )
}
