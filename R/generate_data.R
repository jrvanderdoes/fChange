
#' Title
#'
#' @param resolution
#' @param N
#'
#' @return
#' @export
#'
#' @examples
generateFunctionalIID <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = TRUE)

  t(fiid)
}
