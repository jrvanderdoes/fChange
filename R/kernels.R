#' Bartlett Kernel
#'
#' Kernel where \eqn{max(0,1-|m/h|), h\neq 0} and \eqn{1, h=0}.
#'
#' Note, this function is vectorized.
#'
#' @param x A numeric value at which to evaluate kernel.
#'  It often indicates current lag divided by window.
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @examples
#' bartlett_kernel(-20:20/15)
bartlett_kernel <- function(x) {
  pmax(0, ifelse(is.nan(1 - abs(x)),1,1 - abs(x)) )
}


#' Flat Kernel
#'
#' Kernel where \eqn{min(1,max(1.1-|m/h|,0))}
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references
#' L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional time series", Journal of
#'  Econometrics, 179(1), 66-82.
#'
#' @examples
#' flat_kernel(-20:20/15)
flat_kernel = function(x) {
  pmin(1, pmax(1.1 - abs(x), 0))
}

#' Parzen Kernel
#'
#' Kernel where \eqn{1 - 6 * x^2 + 6 * |x|^3, |x|<=0.5},
#'  \eqn{ 2 * (1 - |x|)^3, 0.5<|x|<1}, and \eqn{0, |x|>1}.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @examples
#' parzen_kernel(-20:20/15)
parzen_kernel <- function(x) {
  ifelse(abs(x) <= 1,
         ifelse(abs(x) <= 0.5,
                1 - 6 * x^2 + 6 * abs(x)^3,
                2 * (1 - abs(x))^3 ),
         0 )
}


#' Daniell Kernel
#'
#' Kernel where \eqn{sin(pi * x) / (pi * x), x!=0} and \eqn{1,X\x=0}.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @examples
#' daniell_kernel(-20:20/15)
daniell_kernel <- function(x) {
  ifelse( x != 0, sin(pi * x) / (pi * x), 1 )
}
