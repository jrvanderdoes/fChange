
#' Generate Data via Ornstein-Uhlenbeck Process
#'
#' Generate autoregressive process with errors according the Ornstein-Uhlenbeck
#'  process.
#'
#' @param N Numeric for the number of observations.
#' @param v Numeric for resolution of data or a vector specifying the
#'  observation points.
#' @param rho Numeric which indicates the dependence on the previous
#'    curve.
#'
#' @return A dfts object for the generated data.
#' @export
#'
#' @examples
#' generate_ornstein_uhlenbeck(N=100,v=20)
generate_ornstein_uhlenbeck <- function(N, v, rho=0){
  if(is.null(rho))
    rho <- 0
  if(length(v)==1)
    v <- seq(0,1,length.out=v)
  r <- length(v)

  # Covariance structure (OU process Cov)
  comat <- matrix(NA, r, r)
  for (i in 1:r){
    comat[i,] <- exp(-abs(v[i]-v)/2)
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = rep(0,r),
                        Sigma = comat,
                        empirical = TRUE)

  data <- data.frame(matrix(0,ncol = nrow(fiid),nrow=ncol(fiid)))
  data[,1] <- fiid[1,]
  for(i in 2:ncol(data)){
    data[,i] <- rho*data[,i-1] + fiid[i,]
  }

  dfts(data)
}
