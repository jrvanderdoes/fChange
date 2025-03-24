#' Generate FAR(1) Data
#'
#' Function to generate data according to FAR(1) process.
#'
#' @param N Numeric for the number of observations.
#' @param resolution Numeric for resolution of data or a vector specifying the
#'  observation points.
#' @param sd Numeric for standard deviation with Brownian motion.
#' @param dependence Numeric which indicates the dependence on the previous
#'    curve.
#' @param drop_first Booolean if first values should be dropped so the data varies
#'  at the first rather than starting at 0 (given that is the observed first point).
#'  Note this will affect the resolution observed.
#'
#' @return dfts object of the data.
#' @export
#'
#' @examples
#' res <- generate_far1(20,10)
generate_far1 <- function(N, resolution, sd=1, dependence=1/2, drop_first=FALSE){
  if(is.null(drop_first))
    drop_first <- TRUE
  if(is.null(dependence))
    dependence <- 0

  # Setup resolution
  if(length(resolution)==1){
    resolution <- seq(0,1,length.out=resolution)
  }

  K <- function(s,t){exp(-(s^2+t^2)/2)}

  # c <- dependence * 1 / sqrt(pracma::quad2d(function(x,y,K){ K(x,y)^2 },
  #                                  xa = 0,xb = 1,ya = 0, yb = 1,K=K) )
  K_mat <- sapply(seq(0,1,length.out=500),function(x,ys,K){K(x,ys)^2}, K=K, y=seq(0,1,length.out=500))
  c <- dependence * 1 / sqrt(dot_integrate(dot_integrate_col(K_mat)))

  w <- generate_brownian_motion(N = N, v = resolution, sd = sd)$data
  if(drop_first){
    w <- w[-1,]
    r <- length(resolution)-1
    res <- resolution[-1]
  }else{
    r <- length(resolution)
    res <- resolution
  }

  X <- data.frame(matrix(0, ncol=N,nrow=r))
  X[,1] <- w[,1]

  for(i in 2:N){
    for(j in 1:r){
      X[j,i] <- c * stats::integrate(function(s,t,X_lag){
        K(s,t)*X_lag
      },lower = 0,upper = 1,t=res[j], X_lag=X[j,i-1])$value + w[j,i]
    }
  }

  dfts(X, fparam = res)
}
