
#' Generate Data via OU Process
#'
#' @param resolution Numeric for data resolution
#' @param N Numeric for data length
#' @param rho Numeric for amount of dependence
#'
#' @return funts object of generated OU data
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' generate_ou(20,100)
.generate_ou <- function(resolution, N, rho=0){
  if(is.null(rho))
    rho <- 0
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

  funts(data)
}
