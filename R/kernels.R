#' Bartlett Kernel
#'
#' A very common Kernel. This is an obvious option in the Welsh-approximation.
#'
#' @param l A numeric value indicating current lag
#' @param h A numeric value indicating smallest lag with no influence
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @examples
#' bartlett_kernel(1, 2)
#' bartlett_kernel(0, 0)
bartlett_kernel <- function(l, h) {
  pmax(0,ifelse(is.nan(1 - abs(l) / h),1,1 - abs(l) / h))
}


#' Title
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
#' flat_kernel(1, 2)
#' flat_kernel(0, 0)
flat_kernel = function(l,h) {
  pmin(1, pmax(1.1 - abs(x/h), 0))
}
