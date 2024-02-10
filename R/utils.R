#' #' Approximate Integral
#' #'
#' #' This (internal) function computes the estimates the integral for
#' #'  evenly spaced observations.
#' #'
#' #' @param y Numeric vector of data
#' #' @param type (Optional) String indicating type of integration. Currently the
#' #'  options are "Rectangle" and "Trapezoidal", both done using Riemann sums.
#' #'  Default is "Trapezoidal".
#' #'
#' #' @return Numeric indicating the estimated integral of the curve
#' #'
#' #' @noRd
#' #'
#' #' @examples
#' #' .approx_int(rep(1, 10), type = "Rectangle")
#' #' .approx_int(seq(0, 1, length.out = 20))
#' #' .approx_int(seq(0, 1, length.out = 20)^2)
#' .approx_int <- function(y, type = "Trapezoidal") {
#'   if (type == "Rectangle") { # Rect Approx
#'     value <- sum(y) / length(y)
#'   } else if (type == "Trapezoidal") { # Trap Approx
#'     value <- sum(y[-1] + y[-length(y)]) / (2 * (length(y) - 1))
#'   } else {
#'     stop("Error: This type is not yet implemented. See documentation.")
#'   }
#'
#'   value
#' }


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


#' Title
#'
#' @param vals
#' @param n
#'
#' @return
#' @export
#'
#' @examples
.select_n <- function(vals,n){
  vals[round(seq(1,length(vals),length.out=n))]
}


#' Get Chunks
#'
#' This (internal) function splits the vector x into a chunksN number of
#'     subsegments. The values of x are kept in order (i.e. no scrambled).
#'
#' @param x Vector of values
#' @param chunksN Numeric indicating the number of chunks to split X into
#'
#' @return A list with chunksN items, each containing an similiar sized subset
#'     of the original vector
#'
#' @noRd
#'
#' @examples
#' .getChunks(1:100, 1)
#' .getChunks(1:100, 2)
#' .getChunks(1:100, 5)
.getChunks <- function(x, chunksN) {
  if (chunksN < 2) {
    return(x)
  }
  split(x, cut(x, chunksN, labels = FALSE))
}


#' Convert List of Samples into a Data Frame
#'
#' This (internal) function takes a list with differ length data.frames or
#'     vectors and pads them all to make a clean data.frame.
#'
#' @param data_list List of elements to be combined to a data.frame.
#'
#' @return Data.frame of the data in data_list
#'
#' @noRd
#'
#' @examples
#' # This is an internal function and will not be viewable to user. See
#' #     generalized_resampling.
.convertSamplesToDF <- function(data_list) {
  m <- length(data_list)
  maxLen <- 0
  for (ii in 1:m) {
    maxLen <- max(maxLen, length(data_list[[ii]]))
  }

  data_df <- data.frame(matrix(nrow = maxLen, ncol = m))

  for (ii in 1:length(data_list)) {
    data_df[, ii] <- c(
      data_list[[ii]],
      rep(NA, maxLen - length(data_list[[ii]]))
    )
  }

  data_df
}
