
#' Potential trim function
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @return Integer value indicating the amount to trim
#' @export
#'
#' @references Rice, G, Zhang, C. (2021), \emph{Consistency of binary
#'     segmentation for multiple change-point estimation with functional data}
#'     (https://arxiv.org/pdf/2001.00093.pdf)
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
trim_function <-function(data){
  data <- as.data.frame(data)
  n <- ncol(data)

  squ_diff <- rep(NA,n-1)
  for(i in 2:n){
    squ_diff[i-1] <- sum((data[,i] - data[,i-1])^2)
  }
  sigma_hat_n_2 <- median(squ_diff)/2

  xi_n <- sqrt(sigma_hat_n_2) * sqrt(3 * log(n))

  round(xi_n)
}
