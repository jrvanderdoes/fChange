#' Covariance Kernel Change
#'
#' This method implements a method for detection of changes in the covariance
#'  kernel of functional data based on the norms of a generally weighted
#'  process of partial sample estimates.
#'
#' Upcoming: Size estimates may be included (already coded).
#'
#' @param X Numeric data.frame of functional data observations--rows for
#'    evaluated values and columns indicating FD
#' @param kappa (Optional) Numeric used for weighting and such. Default is 1/4.
#' @param len (Optional) Numeric for window/repetitions for covariance change.
#'    Default is 30.
#'
#' @return Location of change. NA for no change and numeric if there is a change.
#'
#' @keywords internal
#' @noRd
#'
#' @references Horvath, L., Rice, G., & Zhao, Y. (2022). Change point analysis
#'  of covariance functions: A weighted cumulative sum approach. Journal of
#'  Multivariate Analysis, 189, 104877-.
#'
#' @examples
#' #result <- .change_covariance_kernel(electricity$data[,1:18], len=20)
.change_covariance_kernel <- function(X, statistic, critical,
                                      kappa = 1 / 4, len = 30,
                                      blocksize=1, M=1000, resample_blocks='separate',
                                      replace=TRUE, K=bartlett_kernel) {

  if(statistic=='Tn'){
    tmp <- .covariance_statistic(X$data, statistic=statistic, kappa = kappa, location=TRUE)
    stat <- tmp[1]
    location <- tmp[2]

    if(critical == 'simulation'){
      simulations <- .covariance_simulations(xf = X$data, len = len, kappa = kappa,
                                             M=M, statistic=statistic, K=K)
    } else if(critical == 'resample'){
      simulations <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                           type = resample_blocks, replace = replace,
                           fn = .covariance_statistic,
                           statistic=statistic, kappa = kappa)
    }
  } else{
    stop('Only Tn statistic currently implemented',call. = FALSE)
  }

  return(list(
    'pvalue' = sum(stat <= simulations) / length(simulations),
    'location' = location#,
    # 'statistic' = stat,
    # 'simulations' = simulations
  ))
  # if (stat_d0[[1]] > cv_d0[2]) {
  #   return(stat_d0[[2]])
  #   # kstar = stat_d0[[2]]
  #   # changetau = tau_est(dh1, kstar, len=20)
  #   # cbar = changetau[[1]]
  #   # tau = changetau[[2]]
  #   # changestat=l2norm(cbar)^2/tau*((kstar/samplesize)-truek)
  #   # stat_vec[i]=changestat
  # }
  #
  # return(NA)
}


#' Compute critical values for .covariance_statistic ( \eqn{T_N(\kappa)} )
#'
#' This (internal) function computes the critical values for weighted Tn statistic
#'
#' @inheritParams .covariance_statistic
#' @param len Numeric for window/repetitions for covariance change.
#'
#' @return Numeric critical values (0.9, 0.95, 0.99)
#'
#' @keywords internal
#' @noRd
.covariance_simulations <- function(xf, len, kappa, M, statistic, K=bartlett_kernel) {
  grid_point <- nrow(xf)
  N <- ncol(xf)

  ## cov weight function
  times <- 1:grid_point / grid_point
  wmat <- matrix(NA, grid_point - 2, grid_point - 2)
  for (i in 2:(grid_point - 1)) {
    for (j in 2:(grid_point - 1)) {
      wmat[i - 1, j - 1] <- (min(times[i], times[j]) - times[i] * times[j]) /
        ((times[i] * (1 - times[i]))^kappa * (times[j] * (1 - times[j]))^kappa)
    }
  }
  weig <- as.vector(svd(wmat / grid_point)$d)

  ## cov operators
  rref <- stats::runif(len, 0, 1)
  rref <- c(sort(rref), 1)
  rrefind <- round(rref * dim(xf)[1])
  rrefind[which(rrefind == 0)] <- 1
  xfMC <- xf[rrefind, ]

  # xdm <- apply(xfMC, 2, function(x, xmean) {
  #   x - xmean
  # }, xmean = rowMeans(xfMC))
  xdm <- center(xfMC)

  # zi <- zm <- array(0, c((len + 1), (len + 1), N))
  # for (i in 1:N) {
  #   zi[, , i] <- xdm[, i] %o% xdm[, i]
  # }
  zi <- array(apply(xdm, 2, tcrossprod),
              dim=c(nrow(xdm),nrow(xdm),ncol(xdm)))
  # zimean <-
  #   apply(zi, c(1, 2), mean)
  zimean <- rowMeans(zi,dims = 2)

  zm <- array(0, c((len + 1), (len + 1), N))
  for (i in 1:N) {
    zm[, , i] <- zi[, , i] - zimean
  }

  lrcov <- long_run_covariance_4tensor(zm, K=K)
  lrcov <- tensorA::as.tensor(round(lrcov / (len + 1)^2, 6) )
  eigvals <- tensorA::svd.tensor(lrcov, c(3, 4), by = "e")
  eigmat <- as.vector(eigvals$d)

  lim_sum <- 0
  for (ell in 1:length(eigmat)) {
    klim <- 0
    for (k in 1:length(weig)) {
      Nm <- stats::rnorm(M, mean = 0, sd = 1)
      klim <- klim + eigmat[ell] * weig[k] * Nm^2
    }
    lim_sum <- lim_sum + klim
  }

  lim_sum
}


#' Compute Zn Statistic for Functional Covariance Changes Under Change
#'
#' @inheritParams .ZNstat
#'
#' @return Numeric Zn statistic
#'
#' @keywords internal
#' @noRd
.covariance_statistic_cp <- function(xdm, u) {
  grid_point <- nrow(xdm)
  N <- ncol(xdm)
  k <- floor(N * u)

  prek <- matrix(rowSums(apply(as.matrix(xdm[, 1:k]), 2, function(x) {
    x %o% x
  })), grid_point, grid_point)
  postk <- matrix(rowSums(apply(as.matrix(xdm[, (k+1):N]), 2, function(x) {
    x %o% x
  })), grid_point, grid_point)
  ZNu <- (prek - (k / N) * (prek+postk))

  ZNu
}


#' Compute the Weighted Tn Statistic
#'
#' The function computes the weighted Tn statistic, introduced after Theorem 2.3
#'
#' @param xf Data.frame of numerics. Finite realization of functional time
#'    series data, where curves are stored in columns.
#' @param kappa Numeric used for weighting and such
#'
#' @return List of two values: (stat) giving the weighted Tn statistic and
#'    (changepoint) numeric for change location
#'
#' @keywords internal
#' @noRd
.covariance_statistic <- function(xf, statistic, kappa, location=FALSE) {

  if(statistic !='Tn') stop('Statistic must be `Tn`',call. = FALSE)

  grid_point <- nrow(xf)
  N <- ncol(xf)

  xdm <- apply(xf, 2, function(x, xmean) {
    x - xmean
  }, xmean = rowMeans(xf))

  uind <- seq(0, 1, length = N + 1)[2:(N + 1)]
  zn2 <- list()
  zn_cp <- c(rep(0, N))
  for (i in 1:(N - 1)) {
    Zn_stat <- .covariance_statistic_cp(xdm, uind[i])
    zn2[[i]] <- (Zn_stat/sqrt(N))^2 / ( (uind[i] * (1 - uind[i]))^(2 * kappa) )

    zn_cp[i] <- (N / (i * (N - i)))^kappa *
      dot_integrate(dot_integrate_col( Zn_stat^2))
    #   .int_approx_tensor( (.covariance_statistic_cp(xdm, uind[i]))^2 )
    # sum( sum( (.covariance_statistic_cp(xdm, uind[i]))^2 / grid_point ) / grid_point )
  }
  inm <- Reduce(`+`, zn2) / N
  stat <- (1 / grid_point)^2 * sum(inm)

  if(!location) return(stat)

  # TODO:: Remove Trim
  mcp <- max(zn_cp[(0.1 * N):(0.9 * N)])
  changepoint <- which(zn_cp == mcp)

  c(stat, changepoint)
}


#' Long-run covariance estimator using tensor operator
#'
#' @param dat An array with dimension (grid_point,grid_point,N)
#'
#' @return Numeric matrix for long-run covariance
#'
#' @keywords internal
#' @noRd
long_run_covariance_4tensor <- function(dat, K=bartlett_kernel) {
  grid_point <- dim(dat)[1]
  Tval <- dim(dat)[3]
  datmean <- apply(dat, c(1, 2), mean)
  center_dat <- sweep(dat, 1:2, datmean)

  .cov_l <- function(band, nval, K=bartlett_kernel) {
    cov_sum <- .gamma_l(0, nval)

    for (ik in 1:(nval - 1)) {
      cov_sum <- cov_sum + K(ik / band) * (2 * .gamma_l(ik, nval))
    }

    cov_sum
  }

  .gamma_l <- function(lag, Tval) {
    gamma_lag_sum <- 0
    if (lag >= 0) {
      for (ij in 1:(Tval - lag)) {
        gamma_lag_sum <- gamma_lag_sum + center_dat[, , ij] %o% center_dat[, , (ij + lag)]
      }
    } else {
      for (ij in 1:(Tval + lag)) {
        gamma_lag_sum <- gamma_lag_sum + center_dat[, , (ij - lag)] %o% center_dat[, ij]
      }
    }

    gamma_lag_sum / (Tval - lag)
  }

  hat_h_opt <- Tval^(1 / 4)
  lr_covop <- .cov_l(band = hat_h_opt, nval = Tval, K=K)

  lr_covop
}



# #' L2 Norm
# #'
# #' This (internal) function computes the L2 norm of the data.
# #'
# #' See use in tau_est
# #'
# #' @param vec Vector of numerics to compute L2 norm
# #'
# #' @return Numeric L2 norm value
# #'
# #' @keywords internal
# #' @noRd
# .l2norm <- function(vec) {
#   return(sqrt(sum(vec^2)))
# }

# #' Compute Zn Statistic for Functional Covariance Changes Under Null
# #'
# #' This (internal) function implements Theorem 2.1 / 2.2.
# #'
# #' @param xdm Data.frame of numerics. It is the finite realization of DEMEAN-ed
# #'  functional time series data, where curves are stored in columns.
# #' @param u Numeric in \eqn{(0, 1)}. That is, a fraction index over the interval
# #' \eqn{[0, 1]}.
# #'
# #' @return Numeric Zn statistic
# #'
# #' @keywords internal
# #' @noRd
# .ZNstat <- function(xdm, u) {
#   grid_point <- nrow(xdm)
#   N <- ncol(xdm)
#   k <- floor(N * u)
#
#   prek <- matrix(rowSums(apply(
#     as.matrix(xdm[, 1:k]), 2,
#     function(x) {
#       x %o% x
#     }
#   )), grid_point, grid_point)
#   fullk <- matrix(rowSums(apply(
#     as.matrix(xdm), 2,
#     function(x) {
#       x %o% x
#     }
#   )), grid_point, grid_point)
#   ZNu <- N^(-1 / 2) * (prek - (k / N) * fullk)
#
#   ZNu
# }



# #' Kernels for applying weights
# #'
# #' This (internal) function computes the kernels and related weights
# #'
# #' @param x <Add Information>
# #' @param kernel String indicating kernel to use.
# #' @param normalize (Optional) Boolean indicating if x should be normalized.
# #'    Default is FALSE.
# #'
# #' @return Numeric for k weight
# #'
# #' @keywords internal
# #' @noRd
# .kweights <- function(x,
#                      kernel = c(
#                        "Truncated", "Bartlett", "Parzen", "Tukey-Hanning",
#                        "Quadratic Spectral"
#                      ),
#                      normalize = FALSE) {
#   kernel <- match.arg(kernel)
#   if (normalize) {
#     ca <- switch(kernel,
#       Truncated = 2,
#       Bartlett = 2 / 3,
#       Parzen = 0.539285,
#       `Tukey-Hanning` = 3 / 4,
#       `Quadratic Spectral` = 1
#    )
#   } else {
#     ca <- 1
#   }
#   switch(kernel,
#     Truncated = {
#       ifelse(ca * x > 1, 0, 1)
#     },
#     Bartlett = {
#       ifelse(ca * x > 1, 0, 1 - abs(ca * x))
#     },
#     Parzen = {
#       ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5, 1 - 6 * (ca *
#         x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
#     },
#     `Tukey-Hanning` = {
#       ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x)) / 2)
#     },
#     `Quadratic Spectral` = {
#       y <- 6 * pi * x / 5
#       ifelse(x < 1e-04, 1, 3 * (1 / y)^2 * (sin(y) / y - cos(y)))
#     }
#   )
# }


# #' Approximate Integral of Tensor
# #'
# #' This (internal) function using a Riemann sum to approximate the integral
# #'
# #' @param x  4-dimensional tensor
# #'
# #' @return Approximate integral value
# #'
# #' @keywords internal
# #' @noRd
# .int_approx_tensor <- function(x) {
#   dt <- length(dim(x))
#   temp_n <- nrow(x)
#
#   return(sum(x) / (temp_n^dt))
# }

