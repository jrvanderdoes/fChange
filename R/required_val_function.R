#' Potential required value function
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @return Numeric value indicating necessary value
#' @export
#'
#' @references Rice, G, Zhang, C. (2021), \emph{Consistency of binary
#'     segmentation for multiple change-point estimation with functional data}
#'     (https://www.sciencedirect.com/science/article/pii/S0167715221001905)
#'
#' @examples
#' required_val_function(electricity)
required_val_function <- function(data) {
  data <- as.data.frame(data)
  n <- ncol(data)

  squ_diff <- rep(NA, n - 1)
  for (i in 2:n) {
    squ_diff[i - 1] <- sum((data[, i] - data[, i - 1])^2)
  }
  sigma_hat_n_2 <- stats::median(squ_diff) / 2

  xi_n <- sqrt(sigma_hat_n_2) * sqrt(3 * log(n))

  xi_n
}
