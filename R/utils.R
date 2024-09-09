
#' Check data for functions
#'
#' @param X funts or object that can be directly converted into funts object
#'
#' @return funts object or error if the data is incorrect
#'
#' @examples
#' #.check_data(funts(electricity))
#' #.check_data(electricity)
.check_data <- function(X){
  if(class(X)[1]=='funts') return(X)

  if(class(X)[1]=='matrix' || class(X)[1]=='data.frame') return(funts(X))

  stop('Check format of input data',call. = FALSE)
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


#' Discretization for curve to n observations
#'
#' This will duplicate values if \code{length(vals)>n}.
#'
#' @param vals Vector of values
#' @param n Numeric for number of final observed points / final discretization
#'
#' @return Vector with the data selected at the given level n
#'
#' @noRd
.select_n <- function(vals,n){
  vals[round(seq(1,length(vals),length.out=n))]
}


#' Get Chunks
#'
#' This (internal) function splits the vector x into a \code{chunksN} number of
#'     subsegments. The values of x are kept in order (i.e. not scrambled).
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
  chunksN <- round(chunksN)

  if (chunksN < 2) {
    return(x)
  }
  split(x, cut(x, chunksN, labels = FALSE))
}


#' Convert List of Samples into a Data Frame
#'
#' This (internal) function takes a list with differ length data.frames or
#'  vectors and pads them all to make a clean data.frame.
#'
#' This is an internal function and will not be viewable to user. See
#'  generalized_resampling for usage
#' @param data_list List of elements to be combined to a data.frame.
#'
#' @return Data.frame of the data in data_list
#'
#' @noRd
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
