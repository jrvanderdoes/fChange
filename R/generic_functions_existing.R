

# #' @rdname median
# #'
# #' @param x data
# #' @param na.rm boolean
# #'
# #' @export
# median <- function(x, na.rm = FALSE, ...) UseMethod("median")
# #' @rdname acf
# #'
# #' @export
# median.default <- function(x, na.rm = FALSE, ...) stats::acf(x, na.rm=na.rm, ...)
#' Median for funts
#'
#' @inheritParams max.funts
#'
#' @return Numeric for the median at each resolution
#' @export
#'
#' @importFrom stats median
#'
#' @examples
#' # median(funts(electricity))
median.funts <- function(x, na.rm = FALSE, ...) {
  apply(x$data, MARGIN = 1, stats::median, na.rm=na.rm, ...)
}


# #' @rdname quantile
# #'
# #' @param x data
# #'
# #' @export
# quantile <- function(x, ...) UseMethod("quantile")
# #' @rdname quantile
# #'
# #' @export
# quantile.default <- function(x, ...) stats::quantile(x,...)

#' Quantile funts
#'
#' @param x A funts object
#' @param ... Additional parameters to pass into quantile function
#'
#' @return Matrix with columns for each requested quantile
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' # quantile(funts(electricity))
#' # quantile(funts(electricity),probs = 0.95)
quantile.funts <- function(x, probs = seq(0, 1, 0.25), ...){
  if(length(probs)>1){
    output <- t(apply(x$data, MARGIN = 1, stats::quantile, probs=probs, ...))
  } else{
    output <- data.frame(apply(x$data, MARGIN = 1, stats::quantile, probs=probs, ...))
    colnames(output) <- paste0(probs*100,'%')
  }

  output
}


# #' @rdname lag
# #'
# #' @param x data
# #'
# #' @export
# lag <- function(x, ...) UseMethod("lag")
#' Lag funts
#'
#' @param x funts object
#' @param k integer indicating the number of lags for the data
#' @param ... Unused additional parameters
#'
#' @return A funts object
#' @export
#'
#' @importFrom stats lag
#'
#' @examples
#' # lag.funts(funts(electricity))
lag.funts <- function(x, k, ...) {
  stopifnot(k >= 0)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")
  #
  dat1 <- x$data[,-c(1:round(k))]
  funts(dat1,labels = x$labels[-c(1:round(k))])
}
