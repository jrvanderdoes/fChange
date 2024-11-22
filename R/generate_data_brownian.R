#' Generate a Brownian Motion Process
#'
#' Generate a functional time series according to an iid Brownian Motion process.
#'  Each observation is discretized in the points indicated in \code{v}.
#'
#' @param N Numeric. The number of observations for the generated data.
#' @param v Numeric. Discretization points of the curves. By default it is
#'  evenly spaced on [0,1], i.e. \code{seq(from = 0, to = 1, length.out = 100)}.
#' @param sd Numeric. Standard deviation of the Brownian Motion process.
#'  The default is \code{1}.
#'
#' @return Functional time series object
#' @export
#'
#' @examples
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 20)
#' bmotion <- generate_brownian_motion(N, v, sd = 1)
generate_brownian_motion <- function(
    N, v = seq(from = 0, to = 1, length.out = 20), sd = 1){

  res <- length(v)

  if(N>1){
    data <- apply(t(sapply(c(0,diff(v)),
                           function(d,N,sd){stats::rnorm(N,sd=sd*sqrt(d))},N=N,sd=sd)),
                  MARGIN = 2, cumsum)
  }else{
    data <- cumsum(sapply(c(0,diff(v)),
                           function(d,N,sd){stats::rnorm(N,sd=sd*sqrt(d))},N=N,sd=sd)
                   )
  }

  funts(X=data,intraobs = v)
}


#' Generate a Brownian Bridge Process
#'
#' Generate a functional time series from an iid Brownian Bridge process.
#'  If \eqn{W(t)} is a Wiener process, the Brownian Bridge is
#'  defined as \eqn{W(t) - tW(1)}.
#'  Each functional observation is discretized in the points
#'  indicated in \code{v}.
#'
#' @inheritParams generate_brownian_motion
#'
#' @return Functional time series object
#' @export
#'
#' @examples
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 20)
#' bbridge <- generate_brownian_bridge(N, v, sd=1)
generate_brownian_bridge <- function(
    N, v = seq(from = 0, to = 1, length.out = 100), sd = 1){
  ## TODO:: Convert to other BB (below) for large boost in speed

  res <- length(v)

  # Initialize
  data <- generate_brownian_motion(N, v = v, sd)$data
  data <- data -
    t(data[res,] * t(matrix(rep(v, times = N)/max(v), ncol = N, nrow = res)))

  funts(X=data,intraobs = v)
}


generate_brownian_bridge1 <- function(
    N, v = seq(from = 0, to = 1, length.out = 100), sd = 1){

  n <- length(v)
  dt <- diff(v)
  t <- seq(v[1], v[n], length = n + 1)
  BB <- sapply(1:N, function(m,v,n){
    X <- c(0, cumsum(stats::rnorm(n-1) * sqrt(dt) * sd))
    v[1] + X - (v - v[1])/(v[n] - v[1]) * (X[n] - v[1] + v[1])
  },v=v,n=n)

  funts(X = BB, intraobs = v)
}
