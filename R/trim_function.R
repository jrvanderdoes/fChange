
#' Basic trim function
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param minVal Integer indicating the minimal trim amount
#' @param ... Unused to allow for use in other functions
#'
#' @return Integer value indicating the amount to trim
#' @export
#'
#' @references Rice, G, Zhang, C. (2021), \emph{Consistency of binary
#'     segmentation for multiple change-point estimation with functional data}
#'     (https://www.sciencedirect.com/science/article/pii/S0167715221001905)
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(12,12,12),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,0,1),
#'     distsArray = c('Normal','Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0,0))
#'
#' trim_function(data_KL)
trim_function <- function(data, minVal = 10, ...){
  max(round(minVal), floor(log(ncol(as.data.frame(data)))),na.rm=T)
}
