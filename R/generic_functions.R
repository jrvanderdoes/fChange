#' Generic ACF/PACF
#'
#' @param object
#' @param ...
#'
#' @name acf
#'
#' @return
#' @export
#'
#' @examples
#' acf(1:10)
NULL


#' @rdname acf
#' @export
acf <- function(object, ...) UseMethod("acf")
#' @rdname acf
#' @export
acf.default <- function(object, ...) stats::acf(object)


#' @rdname acf
#' @export
pacf <- function(object, ...) UseMethod("pacf")
#' @rdname acf
#' @export
pacf.default <- function(object, ...) stats::pacf(object)


#' ACF for Functional Data
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' acf(funts(electricity))
acf.funts <- function(x, ...){
  invisible(.compute_FACF(x$data, x$intraobs, ...))
  # lag.max = NULL, ci=0.95, mc_estimation = TRUE,figure = TRUE, ...)
}


#' PACF for Functional Data
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' pacf(funts(electricity))
pacf.funts <- function(x, ...){
  invisible(.compute_FPACF(x$data, x$intraobs, ...))
    # n_harm = NULL, lag.max = NULL, ci=0.95, figure = TRUE, ...)
}
