#' Generate FAR(1) Data
#'
#' Function to generate date according to FAR(1) process.
#'
#' @param resolution Numeric for resolution of data
#' @param N Numeric for the number of observations
#' @param d Numeric (0,Inf) which indicates the dependence on the previous
#'    curve.
#' @param burnin Numeric for burn-in before data is computed
#'
#' @return Data.frame for functional observations
#' @export
#'
#' @examples
#' res <- generate_far1(24,200)
generate_far1 <- function(resolution, N, d=1/2, burnin=1000){

  times <- seq(0,1,length.out=resolution)
  K <- function(s,t){exp(-(s^2+t^2)/2)}

  c <- d * (1 / sqrt(cubature::adaptIntegrate(
    function(x){ K(x[1],x[2])^2 },
    c(0, 0), c(1, 1))$integral))

  X <- data.frame(matrix(0, ncol=N+burnin,nrow=resolution))
  w  <- sapply(1:(N+burnin), function(x, resolution){
    as.numeric(sde::BM(N=resolution))[-1]
  }, resolution=resolution)
  X[,1] <- w[,1]

  for(i in 2:(N+burnin)){
    for(j in 1:length(times)){
      X[j,i] <- c * stats::integrate(function(s,t,X_lag){
        K(s,t)*X_lag
      },lower = 0,upper = 1,t=times[j], X_lag=X[j,i-1])$value + w[j,i]
    }
  }

  funts(X[,burnin + 1:N])
}
