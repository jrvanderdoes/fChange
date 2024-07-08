#' Basic trim function
#'
#' Basic trim function to select how much data to consider.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param minVal Integer indicating the minimal trim amount
#' @param ... Unused, but exist to allow for use in other functions
#'
#' @return Integer value indicating the amount to trim
#' @export
#'
#' @references Rice, G, Zhang, C. (2021), \emph{Consistency of binary
#'     segmentation for multiple change-point estimation with functional data}
#'     (https://www.sciencedirect.com/science/article/pii/S0167715221001905)
#'
#' @examples
#' trim_function(electricity)
trim_function <- function(data, minVal = 10, ...) {
  max(round(minVal), floor(log(ncol(as.data.frame(data)))), na.rm = TRUE)
}
