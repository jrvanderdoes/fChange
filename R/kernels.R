
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
#' bartlett_kernel(1,2)
#'
#' # Use in Welsh Approximation
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,1),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' welsh_approximation(data_KL,K=bartlett_kernel)
bartlett_kernel <- function(l, h){
  val <- 0
  if(h==l && h==0){
    val <- 1
  }else if(abs(l) < h){
    val <- (1 - abs(l) / h)
  }

  val
}
