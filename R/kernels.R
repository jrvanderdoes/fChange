
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
#' bartlett_kernel(0,0)
bartlett_kernel <- function(l, h){
  val <- 0
  if(h==l && h==0){ # if both 0
    val <- 1
  }else if(abs(l) < h){ # <= no diff as 0
    val <- (1 - abs(l) / h)
  }

  val
}
