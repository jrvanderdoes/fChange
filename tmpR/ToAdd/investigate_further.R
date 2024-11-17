#' Variance of Measurement Error
#'
#' @param X funts obejct
#'
#' @return
#' @export
#'
#' @references Wang, J.-L., Chiou, J.-M., & Mueller, H.-G. (2015). Review
#'  of Functional Data Analysis. \url{https://doi.org/10.48550/arxiv.1507.05135}
#'
#' @examples
#' var_measurement_error(electricity)
var_measurement_error <- function(X){
  # TODO:: Read Functional Data Analysis for Sparse Longitudinal Data
  X <- .check_data(X)

  1/ncol(X) * ( rowSums( (X$data - mean(X))^2 ) - diag(var(X)) )
}

# Functional Correlation , models.
## Study https://arxiv.org/pdf/1507.05135
