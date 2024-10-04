#' Scalar Characteristic Function Change Point Analysis
#'
#' This implements a scalar characteristic change method.
#'
#' @param Y Scalar data
#' @param gam Numeric for weights
#' @param nSims Numeric for number of simulations
#' @param alpha_val Signficance of test
#'
#' @return XXXXX
#' @export
#'
#' @references Hušková, M., & Meintanis, S.G. (2006). Change Point Analysis
#'  based on Empirical Characteristic Functions. Metrika, 63, 145-168.
#'
#' @examples
#' scalarDetection(c(rnorm(50),rnorm(75,mean = 1)))
#' \donttest{
#' scalarDetection(as.numeric(colMeans(electricity)))
#' }
scalarDetection <- function(Y, gam = 0.5, nSims = 200,
                            alpha_val = 0.05) {
  n <- length(Y)

  Tn_s <- rep(NA, n - 1)
  for (k in 1:(n - 1)) {
    Tn_s[k] <- ((k * (n - k)) / n^2)^gam * ((k * (n - k)) / n) *
      stats::integrate(
        function(t, Y, k) {
          abs(.phi_k(Y, t, k) - .phi_k0(Y, t, k))^2 * .w(t)
        },
        lower = 0, upper = 1, Y = Y, k = k
      )[[1]]
  }
  Tn <- max(Tn_s)

  Tn_permute <- sapply(rep(NA, nSims),
    function(tmp, n, gam, Y) {
      Tn_tmp_s <- sapply(1:(n - 1), function(k, n, gam, Y) {
        ((k * (n - k)) / n^2)^gam * ((k * (n - k)) / n) *
          stats::integrate(
            function(t, Y, k) {
              abs(.phi_k(Y, t, k) - .phi_k0(Y, t, k))^2 * .w(t)
            },
            lower = 0, upper = 1, Y = Y[sample(1:length(Y))], k = k
          )[[1]]
      }, n = n, gam = gam, Y = Y)

      max(Tn_tmp_s)
    },
    n = n, gam = gam, Y = Y
  )

  list(
    "dat" = Y,
    "Tn" = Tn,
    "Tn_permute" = Tn_permute,
    "pval" = 1 - stats::ecdf(Tn_permute)(Tn),
    "Loc" = which.max(Tn_s),
    "loc_if_pval" = ifelse(1 - stats::ecdf(Tn_permute)(Tn) < alpha_val,
      which.max(Tn_s), NA
    )
  )
}


#' Compute Phi_k for scalar detection
#'
#' @inheritParams scalarDetection
#' @param t Time dimension
#' @param k potential change
#'
#' @return Numeric
#'
#' @noRd
.phi_k <- function(Y, t, k) {
  sumVal <- 0
  for (j in 1:k) {
    sumVal <- sumVal +
      exp(complex(real = 0, imaginary = 1) * t * Y[j])
  }

  1 / k * sumVal
}


#' Compute Phi_k0 for scalar detection
#'
#' @inheritParams .phi_k
#'
#' @return Numeric
#'
#' @noRd
.phi_k0 <- function(Y, t, k) {
  n <- length(Y)
  sumVal <- 0

  for (j in (k + 1):n) {
    sumVal <- sumVal +
      exp(complex(real = 0, imaginary = 1) * t * Y[j])
  }

  1 / (n - k) * sumVal
}


#' Weighting Scheme for scalar change point detection
#'
#' @inheritParams .phi_k
#' @param a (Optional) Weight to consider. Default is 1.
#'
#' @return Numeric weight
#'
#' @noRd
.w <- function(t, a = 1) {
  exp(-a * abs(t))
}
