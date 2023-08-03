#' Approximate Integral
#'
#' This (internal) function computes the estimates the integral for
#'  evenly spaced observations.
#'
#' @param y Numeric vector of data
#' @param type (Optional) String indicating type of integration. Currently the
#'  options are "Rectangle" and "Trapezoidal", both done using Riemann sums.
#'  Default is "Trapezoidal".
#'
#' @return Numeric indicating the estimated integral of the curve
#'
#' @noRd
#'
#' @examples
#' .approx_int(rep(1, 10), type = "Rectangle")
#' .approx_int(seq(0, 1, length.out = 20))
#' .approx_int(seq(0, 1, length.out = 20)^2)
.approx_int <- function(y, type = "Trapezoidal") {
  if (type == "Rectangle") { # Rect Approx
    value <- sum(y) / length(y)
  } else if (type == "Trapezoidal") { # Trap Approx
    value <- sum(y[-1] + y[-length(y)]) / (2 * (length(y) - 1))
  } else {
    stop("Error: This type is not yet implemented. See documentation.")
  }

  value
}


#' Specify Decimal
#'
#' This (internal) function returns a string of the numbers with the specified
#'  level of decimals (i.e. will add trailing zeroes as needed).
#'
#' @param x Numeric(s) to specify the decimal for
#' @param k Numeric integer indicating the number of decimals to return for each
#'  numeric in x.
#'
#' @return A vector of strings relating the the values in x, but with the
#'  specified number of decimals.
#'
#' @noRd
.specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall = k))
}
