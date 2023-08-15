#' Estimate Long Run Variance
#'
#' This (internal) function estimates CHat
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param v1 Vector of length nrow(X) indicating a random vector in the space
#' @param v2 Vector of length nrow(X) indicating a random vector in the space
#' @param TVal Numeric value indicate the number of FDs.
#' @param h Numeric value indicating the amount of blocking.
#' @param K Function for the Kernel.
#'
#' @return Numeric indicated estimated CHat value
#'
#' @noRd
#'
#' @examples
#' # This is an internal function, see welch_approximation for usage.
.estimateLRV <- function(X, v1, v2, TVal, h, K) {
  eps1 <- exp(1/nrow(X) * complex(real = 0, imaginary = 1) * (t(X) %*% v1) )
  # eps1Bar <- mean(eps1)#1/TVal * sum(eps1)

  eps2 <- exp(1/nrow(X) * complex(real = 0, imaginary = 1) * (t(X) %*% v2) )
  # eps2Bar <- mean(eps2)#1/TVal * sum(eps2)

  gammaHat <- rep(NA, length((1 - TVal):(TVal - 1)))

  for (l in (1 - TVal):(TVal - 1)) {
    gammaHat[TVal + l] <- .getAutocov(
      eps1 = eps1, eps2 = eps2,
      # eps1Bar=eps1Bar, eps2Bar=eps2Bar,
      l = l, TVal = TVal
    )
  }

  .computeLRVEstimate(h = h, TVal = TVal, gammaHat = gammaHat, K = K)
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
#'
#' @examples
#' # This is an internal function, see .estimateLRV for usage.
.getAutocov <- function(eps1, eps2, l, TVal) { # eps1Bar, eps2Bar, l, TVal){
  gammaVal <- 0
  eps1Bar <- mean(eps1)
  eps2Bar <- mean(eps2)

  if (l >= 0) {
    for (j in (1 + l):TVal) {
      gammaVal <- gammaVal +
        (eps1[j - l] - eps1Bar) * Conj(eps2[j] - eps2Bar)
    }
  } else {
    for (j in (1 - l):TVal) {
      gammaVal <- gammaVal +
        (eps1[j + l] - eps1Bar) * Conj(eps2[j] - eps2Bar) # Confirm this j+l swap
    }
  }

  1 / TVal * gammaVal
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
#'
#' @examples
#' # This is an internal function, see .estimateLRV for usage.
.computeLRVEstimate <- function(h, TVal, gammaHat, K) {
  CVal <- 0

  for (l in (1 - TVal):(TVal - 1)) {
    CVal <- CVal + K(l, h) * gammaHat[TVal + l]
  }

  CVal
}
