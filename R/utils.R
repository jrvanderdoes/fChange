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

#' Bootstrap Data
#'
#' @param X
#' @param blockSize
#' @param M
#' @param type
#' @param replace
#'
#' @return
#'
#' @examples
.bootstrap <- function(X, blockSize, M=1000, type='overlapping', replace=TRUE){
  X <- .check_data(X)
  N <- ncol(X$data)

  # Get groups
  if(type=='seperate'){
    idxGroups <- .getChunks(1:ncol(X$data), ncol(X$data) / blockSize)
  } else if(type=='overlapping'){
    idxGroups <- sapply(0:(N-blockSize), function(x,blockSize){
      x + 1:blockSize }, blockSize=blockSize, simplify = FALSE)
  }

  # Get BS sample indices
  idxs <- sapply(1:M, function(i, m, indxs, replace,size,N) {
    samps <- sample(x=1:m, size=size, replace = replace)
    unlist(indxs[samps], use.names = FALSE)[1:N]
  }, m = length(idxGroups), indxs = idxGroups, replace = replace,
  size=ceiling(ncol(X$data) / blockSize), N=N)

  if (!(is.matrix(idxs) | is.data.frame(idxs))) {
    idxs <- .convertSamplesToDF(idxs)
  }

  ## Sample via bootstrap
  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter, Xdata) {
                        Xdata[, stats::na.omit(loop_iter)]
                      }, Xdata=X$data,simplify = F
  )

  bssamples
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
