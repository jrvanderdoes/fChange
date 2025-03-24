#' Generate a Brownian Motion Process
#'
#' Generate a functional time series according to an iid Brownian Motion process.
#'  Each observation is discretized on the points indicated in \code{v}.
#'
#' @param N Numeric. The number of observations for the generated data.
#' @param v Numeric (vector). Discretization points of the curves.This can be
#'  the evaluated points or the number of evenly spaced points on \[0,1\].
#'  By default it is evenly spaced on \[0,1\] with 30 points.
#' @param sd Numeric. Standard deviation of the Brownian Motion process.
#'  The default is \code{1}.
#'
#' @return Functional time series (dfts) object.
#' @export
#'
#' @examples
#' bmotion <- generate_brownian_motion(N=100,
#'   v=c(0,0.25,0.4,0.7, 1, 1.5), sd = 1)
#' bmotion1 <- generate_brownian_motion(N=100,
#'   v=10, sd = 2)
generate_brownian_motion <- function(N, v = 30, sd = 1){

  # Convert resolution to equally spaced points
  if(length(v)==1){
    v <- seq(0,1,length.out=v)
  }
  res <- length(v)

  # Generate data
  if(N>1){
    data <- apply(t(sapply(c(0,diff(v)),
                           function(d,N,sd){stats::rnorm(N,sd=sd*sqrt(d))},N=N,sd=sd)),
                  MARGIN = 2, cumsum)
  }else{
    data <- cumsum(sapply(c(0,diff(v)),
                           function(d,N,sd){stats::rnorm(N,sd=sd*sqrt(d))},N=N,sd=sd)
                   )
  }

  dfts(X=as.matrix(data),fparam = v)
}


#' Generate a Brownian Bridge Process
#'
#' Generate a functional time series from an iid Brownian Bridge process.
#'  If \eqn{W(t)} is a Wiener process, the Brownian Bridge is
#'  defined as \eqn{W(t) - tW(1)}. Each functional observation is discretized on
#'  the points indicated in \code{v}.
#'
#' @inheritParams generate_brownian_motion
#'
#' @return Functional time series (dfts) object.
#' @export
#'
#' @examples
#' bbridge <- generate_brownian_bridge(N=100, v=c(0,0.2,0.6,1,1.3), sd=2)
#' bbridge <- generate_brownian_bridge(N=100, v=10, sd=1)
generate_brownian_bridge <- function(N, v = 30, sd = 1){

  # Convert resolution to equally spaced points
  if(length(v)==1){
    v <- seq(0,1,length.out=v)
  }
  res <- length(v)

  # Generate Data
  data <- generate_brownian_motion(N, v = v, sd)$data
  data <- data -
    t(data[res,] * t(matrix(rep(v, times = N)/max(v), ncol = N, nrow = res)))

  dfts(X=data,fparam = v)
}


#' Second-Order Brownian Bridge
#'
#' @inheritParams generate_brownian_motion
#'
#' @return dfts object
#'
#' @noRd
#' @keywords internal
.generate_brownian_bridge_second_order <- function(
    N, v = seq(from = 0, to = 1, length.out = 100), sd = 1){

  if(length(v)==1){
    v <- seq(0,1,length.out=v)
  }

  n <- length(v)
  dt <- diff(v)
  t <- seq(v[1], v[n], length = n + 1)

  BB <- sapply(1:N, function(m,v,n){
    X <- c(0, cumsum(stats::rnorm(n-1) * sqrt(dt) * sd))
    v[1] + X - (v - v[1])/(v[n] - v[1]) * (X[n] - v[1] + v[1])
  },v=v,n=n)

  dfts(X = BB, fparam = v)
}


#' Generate Second-Level Brownian Motion
#'
#' Generates a second-order Brownian Motion.
#'
#' @param M Numeric number of simulations
#' @param v Observation points
#'
#' @return dfts object of the Brownian Motion
#'
#' @references MacNeill, I. B. (1978). Properties of Sequences of Partial Sums
#'  of Polynomial Regression Residuals with Applications to Tests for Change of
#'  Regression at Unknown Times. The Annals of Statistics, 6(2), 422-433.
#'
#' @keywords internal
#' @noRd
.generate_second_level_brownian_bridge <- function(N,v, sd=1){
  W <- generate_brownian_motion(N,v=v,sd = sd)
  # W <- generate_data(
  #   general=list('resolution'=v,'burnin'=0),
  #   parameters =list(list('N'=N, 'process'='bmotion', 'sd'=sd))
  # )
  V <- W$data +
    (2*v-3*v^2) %*% t(W$data[nrow(W$data),]) +
    (-6*v + 6*v^2) %*% t(dot_integrate_col(W$data, v))

  dfts(X=V,fparam = v)
}
