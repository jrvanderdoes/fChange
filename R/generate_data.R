
#' Title
#'
#' @param resolution
#' @param N
#'
#' @return
#' @export
#'
#' @examples
generateOU <- function(resolution, N, burnin=500, rho=0){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    #comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    comat[i,] <- exp(-abs(times[i]-times)/2)
  }

  fiid <- MASS::mvrnorm(n = N+burnin,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = TRUE)

  data <- data.frame(matrix(0,ncol = nrow(fiid),nrow=ncol(fiid)))
  data[,1] <- fiid[1,]
  for(i in 2:ncol(data)){
    data[,i] <- rho*data[,i-1] + fiid[i,]
  }

  data
}
