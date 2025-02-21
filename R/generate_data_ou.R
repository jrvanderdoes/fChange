
#' Generate Data via Ornstein-Uhlenbeck Process
#'
#' @param resolution Numeric for data resolution or specific observed intratime
#' @param N Numeric for data length
#' @param rho Numeric for amount of dependence
#'
#' @return dfts object of generated OU data
#' @export
#'
#' @examples
#' generate_ornstein_uhlenbeck(20,100)
generate_ornstein_uhlenbeck <- function(N, resolution, rho=0){
  if(is.null(rho))
    rho <- 0
  if(length(resolution)==1)
    resolution <- seq(0,1,length.out=resolution)
  r <- length(resolution)

  # Covariance structure (OU process Cov)
  comat <- matrix(NA, r, r)
  for (i in 1:r){
    comat[i,] <- exp(-abs(resolution[i]-resolution)/2)
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
