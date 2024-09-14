#' Brownian Motion
#' 
#' Creates at J x N matrix, containing N independent Brownian motion sample paths in
#' each of the columns.
#'
#' @param N the number of independent Brownian motion sample paths to compute.
#' @param J the number of steps observed for each sample path (the resolution of the data).
#' 
#' @return  A J x N matrix containing Brownian motion functional data in the columns.
#' @export
#' 
#' @examples
#' b <- brown_motion(250, 50)
brown_motion <- function(N, J) {
  motion <- matrix(nrow = J, ncol = N)
  for (i in 1:N) {
    motion[,i] <- as.vector(sde::BM(N=J-1))
  }
  
  motion
}


#' Functional GARCH(1,1)
#'
#' Simulates an fGARCH(1,1) process with N independent observations, each observed
#'  discretely at J points on the interval [0,1]. Uses the Ornstein-Uhlenbeck process.
#'
#' @param N the number of fGARCH(1,1) curves to sample.
#' @param J the number of points at which each curve is sampled (the resolution of the data).
#' @param delta a parameter used in the variance recursion of the model.
#' @param burn_in the number of initial samples to burn (discard).
#' 
#' @return A list containing two J x N matrices, the former containing the sample of fGARCH(1,1)
#' curves and the latter containing the respective variance values.
#' @export
#' 
#' @examples
#' f <- fgarch11(100, 50)
fgarch11 <- function(N, J, delta=0.01, burn_in=50) {
  grid <- (1:J) / J
  error_mat <- matrix(nrow=J, ncol=N+burn_in)
  
  covariance_mat <- matrix(0, nrow = J, ncol = J)
  for (i in 1:J) {
    covariance_mat[i,] <- exp(-abs(grid[i]-grid) / 2)
  }
  means <- rep(0,J)
  for (i in 1:(N+burn_in)) {
    error_mat[,i] = MASS::mvrnorm(1, mu=means, Sigma=covariance_mat)
  }
  
  alpha_op <- beta_op <- function(t,s) { 12*t*(1-t)*s*(1-s) }
  garch_mat <- sigma2_mat <- matrix(nrow=J, ncol=N+burn_in)
  int_approx <- function(x) {
    sum(x) / NROW(x)
  }
  sigma2_mat[,1] <- rep(delta, J)
  garch_mat[,1] <- sqrt(delta) * error_mat[,1]
  for (i in 2:(N+burn_in)) {
    for (u in 1:J) {
      alpha_op_vec <- alpha_op(grid[u], grid) * (garch_mat[,i-1] ^ 2)
      beta_op_vec <- beta_op(grid[u], grid) * sigma2_mat[,i-1]
      sigma2_mat[u,i] <- delta + int_approx(alpha_op_vec) + int_approx(beta_op_vec)
    }
    garch_mat[,i] <- sqrt(sigma2_mat[,i]) * error_mat[,i]
  }
  
  garch_mat[,burn_in + 1:N]
}


#' Simulation FAR(1,S)
#' 
#' Simulates an FAR(1,S)-fGARCH(1,1) process with N independent observations, each
#'  observed discretely at J points on the interval [0,1].
#'
#' @param N the number of fGARCH(1,1) curves to sample.
#' @param J the number of points at which each curve is sampled (the resolution of the data).
#' @param S the autoregressive operator of the model, between 0 and 1, indicating the level of
#' conditional heteroscedasticity.
#' @param type the assumed model of the error term. The default argument is 'IID', under which
#' the errors are assumed to be independent and identically distributed. The alternative argument
#' is 'fGARCH', which will assume that the errors follow an fGARCH(1,1) process.
#' @param burn_in the number of initial samples to burn (discard).
#' 
#' @return A J x N matrix containing FAR(1,S) functional data in the columns.
#' @export
#' 
#' @examples
#' f <- far1S(100, 50, 0.75)
far1S <- function(N, J, S, type='IID', burn_in=50) {
  grid <- (1:J) / J
  if (type == 'IID') {
    error_mat <- brown_motion(N + burn_in, J)
  } else if (type == 'fGARCH') {
    error_mat <- fgarch11(N + burn_in, J)
  }
  func <- function(t,s) { exp(-((t^2 + s^2) / 2)) }
  sum1 <- 0
  for (t in grid) {
    sum1 <- sum1 + sum(func(t,grid)^2)
  }
  phi_norm_w_out_c <- sqrt(sum1 / (J^2))
  abs_c <- S / phi_norm_w_out_c
  phi_c_t_s <- function(t, s, c=abs_c) {
    c * exp(-((t^2 + s^2) / 2))
  }
  
  far_mat <- matrix(0, nrow=J, ncol=N+burn_in)
  far_mat[,1] <- error_mat[,1]
  for (i in 2:(N+burn_in)) {
    for (j in 1:J) {
      far_mat[j,i] <- (sum(phi_c_t_s(grid[j], grid) * far_mat[,i-1]) / J) + error_mat[j,i]
    }
  }
  
  far_mat[,burn_in + 1:N]
}
