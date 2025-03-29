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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
.getChunks <- function(x, chunksN) {
  chunksN <- round(chunksN)

  if (chunksN < 2) {
    return(x)
  }
  split(x, cut(x, chunksN, labels = FALSE))
}


#' Bootstrap Data
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param blocksize Numeric indicating size of blocks
#' @param M Numeric indicating the number of iterations
#' @param type String for 'overlapping' or 'separate' block bootstrapping
#' @param replace Boolean if data should be sample with or without replacement
#'
#' @return List of permuted data
#'
#' @keywords internal
#' @noRd
.bootstrap <- function(X, blocksize, M=1000, type='overlapping', replace=TRUE,
                       fn=NULL, ...){
  X <- dfts(X)
  N <- ncol(X$data)

  # Get groups
  if(type=='separate'){
    idxGroups <- .getChunks(1:ncol(X$data), ncol(X$data) / blocksize)
  } else if(type=='overlapping'){
    idxGroups <- sapply(0:(N-blocksize), function(x,blocksize){
      x + 1:blocksize }, blocksize=blocksize, simplify = FALSE)
  }

  # Get BS sample indices
  #   Early code to manage worst case of picking smallest group the entire time
  #   Only can occur with replacement (otherwise will pick all)
  min_len <- min(lengths(idxGroups))
  size <- ceiling(ncol(X$data) / blocksize)
  size_equal_len <- size + ceiling((ncol(X$data)-min_len*size)/min_len)*replace

  idxs <- sapply(1:M, function(i, m, indxs, replace,size, N) {
    samps <- sample(x=1:m, size=size, replace = replace)
    unlist(indxs[samps], use.names = FALSE)[1:N]
  }, m = length(idxGroups), indxs = idxGroups, replace = replace,
  size=size_equal_len, N=N)



  # if (!(is.matrix(idxs) | is.data.frame(idxs))) {
  #   stop('This is an internal error, please report', call. = FALSE)
  #   idxs <- .convertSamplesToDF(idxs)
  # }

  ## Sample via bootstrap and return
  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter, Xdata) {
                        Xdata[, stats::na.omit(loop_iter)]
                      }, Xdata=X$data, simplify = FALSE
  )

  ## Apply fn if given
  if(!is.null(fn)){
    result <- unlist(lapply(bssamples,fn, ...))
  } else{
    result <- bssamples
  }

  result
}


# #' Convert List of Samples into a Data Frame
# #'
# #' This (internal) function takes a list with differ length data.frames or
# #'  vectors and pads them all to make a clean data.frame.
# #'
# #' This is an internal function and will not be viewable to user. See
# #'  generalized_resampling for usage
# #' @param data_list List of elements to be combined to a data.frame.
# #'
# #' @return Data.frame of the data in data_list
# #'
# #' @keywords internal
# #' @noRd
# .convertSamplesToDF <- function(data_list) {
#   m <- length(data_list)
#   maxLen <- 0
#   for (ii in 1:m) {
#     maxLen <- max(maxLen, length(data_list[[ii]]))
#   }
#
#   data_df <- data.frame(matrix(nrow = maxLen, ncol = m))
#
#   for (ii in 1:length(data_list)) {
#     data_df[, ii] <- c(
#       data_list[[ii]],
#       rep(NA, maxLen - length(data_list[[ii]]))
#     )
#   }
#
#   data_df
# }
