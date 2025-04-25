#' Functional Spherical Autocorrelation Function
#'
#' This function offers a graphical summary of the fSACF of a
#'  functional time series (FTS) across different time lags \eqn{h = 1:H}.
#'  It also plots \eqn{100 \times (1-\alpha)\%}  confidence bounds developed
#'  under strong white noise (SWN) assumption for all lags \eqn{h = 1:H}.
#'
#' This function computes and plots functional spherical
#'  autocorrelation coefficients at lag \eqn{h}, for \eqn{h = 1:H}.
#'  The fSACF at lag \eqn{h} is computed by the average of the inner product of
#'  lagged pairs of the series \eqn{X_i} and \eqn{X_{i+h}} that have been
#'  centered and scaled:
#'  \deqn{
#'   \tilde\rho_h=\frac{1}{N}\sum_{i=1}^{N-h} \langle \frac{X_i - \tilde{\mu}}{\|X_i - \tilde{\mu}\|}, \frac{X_{i+h} - \tilde{\mu}}{\|X_{i+h} - \tilde{\mu}\|} \rangle,\ \ \ \ 0 \le h < N,
#'  }
#'  where \eqn{\tilde{\mu}} is the estimated spatial median of the series.
#'  It also computes estimated asymptotic \eqn{(1-\alpha)100 \%} confidence lower
#'  and upper bounds, under the SWN assumption.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param lag.max A positive integer value. The maximum lag for which to compute
#'  the coefficients and confidence bounds.
#' @param alpha Significance in \[0,1\] for intervals when forecasting.
#' @param figure Logical. If \code{TRUE}, prints plot for the estimated
#'  function with the specified bounds.
#'
#' @return List with sACF values and plots
#'
#' @export
#'
#' @seealso [acf()]
#'
#' @references Yeh C.K., Rice G., Dubin J.A. (2023). Functional spherical
#'  autocorrelation: A robust estimate of the autocorrelation of a functional
#'  time series. Electronic Journal of Statistics, 17, 650–687.
#'
#' @examples
#' sacf(electricity)
#' sacf(generate_brownian_motion(100) )
sacf <- function (X, lag.max = 20, alpha = 0.05, figure = TRUE) {
  # Test Input
  if ((lag.max < 1) | (lag.max%%1 != 0)) {
    stop("The parameter 'lag.max' must be a positive integer.")
  }
  if ((alpha > 1) | (alpha < 0)) {
    stop("The 'alpha' parameter must be a value between 0 and 1.")
  }

  X <- dfts(X)
  N <- ncol(X)
  lags <- 1:lag.max


  # Setup Data
  res_raw <- .calculate_sacf(X, H = lag.max)
  # res_raw <- .calculate_sacf(obs = t(f_data), H)

  # coefficients = rep(0, H)
  # B_iid_bounds_L = rep(0, H)
  # B_iid_bounds_U = rep(0, H)

  res_raw <- as.matrix(res_raw)

  # Compute Coeffifients and Bounds
  coefficients <-  res_raw[,3]
  B_iid_bounds <- stats::qt(alpha/2, N - 1) * res_raw[,5]/sqrt(N)
  # B_iid_bounds_U <- stats::qt(1 - alpha/2, N - 1) * res_raw[,5]/sqrt(N)

  # Plotting

  # plot(lags, coefficients, ylim = c(min(coefficients, B_iid_bounds_L[1])-0.05,
  #                                   max(coefficients, B_iid_bounds_U[1])+0.3),
  #      type = 'h', xlab = "Lag", ylab = "Coefficient",
  #      main = "fSACF plot")
  # lines(B_iid_bounds_L, col = 'red', lty = 'solid')
  # lines(B_iid_bounds_U, col = 'red', lty = 'solid')
  # lines(rep(0,H), col = 'black', lty = 'solid')
  #
  # legend('topleft', legend = c('Estimated Spherical Autocorrelation Coefficients',
  #                              paste('SWN ', (1 - alpha) * 100,'% Confidence Bound', sep=' ')),
  #        col = c('black', 'red'),
  #        lty = c("solid", "solid"), cex=0.9, y.intersp = 1, bty = "n")
  #

  plt <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(
      x = lags,
      y= pmin(0,coefficients), yend=pmax(0,coefficients)),
      # position = ggplot2::position_dodge2(preserve='single'),
      # stat = 'identity',
      col = 'black', #width=0.2, fill='darkgray'
    ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0), col="black") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(size=24),
                   axis.text = ggplot2::element_text(size=20)) +
    ggplot2::xlab('Lag') +
    ggplot2::ylab(NULL) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=-B_iid_bounds), col="#0073C2FF",
                        linetype='dashed') +
    ggplot2::geom_hline(ggplot2::aes(yintercept=B_iid_bounds), col="#0073C2FF",
                        linetype='dashed')

  if(figure){
    print(plt)
  }

  invisible( list('acfs'=coefficients, 'SWN'=B_iid_bounds, #'WWN'=WWN_bound,
                  'plot'=plt) )
}


#' Temporary Function for SACF
#'
#' @inheritParams sacf
#' @param H See max.lag
#'
#' @returns Data.frame with info
#'
#' @noRd
#' @keywords internal
.calculate_sacf <- function(X, H){
  N <- ncol(X)#nrow(obs)
  r <- nrow(X)#ncol(obs)

  norm_raw <- .calc_norm_matrix(t(X$data))
  # norm_raw <- .calc_norm_matrix(obs)

  SpMedian <- as.numeric(.spatial_median(seq(0, 1, length.out = r), t(X$data))$med)
  # SpMedian <- as.numeric(.spatial_median(seq(0, 1, length.out = Nt), obs)$med)

  obs_center <- t(
    sapply(1:N, function(i){
      X$data[, i] - SpMedian
    }) )
  # obs_center1 <- t(
  #   sapply(1:N, function(i){
  #     obs[i, ] - SpMedian
  #   }) )

  norm_center <- .calc_norm_matrix(obs_center)
  # norm_center <- .calc_norm_matrix(obs_center1)

  rho_raw <- sapply(1:H, FUN = function(h){
    res <- sapply(1:(N-h), function(i){
      sum(t(X$data)[i, ]/norm_raw[i] * t(X$data)[i+h, ]/norm_raw[i+h])/(r-1)
    })
    sum(res)/N
  })

  rho_cen <- sapply(1:H, FUN = function(h){
    res <- sapply(1:(N-h), function(i){
      sum(obs_center[i, ]/norm_center[i] * obs_center[i+h, ]/norm_center[i+h])/(r-1)
    })
    sum(res)/N
  })

  my_var_raw <- .calc_var(t(X$data)/norm_raw) # .calc_var(obs/norm_raw)
  my_var_cen <- .calc_var_raw(t(X$data)) # .calc_var_raw(obs)

  cbind(1:H, rho_raw,  rho_cen,
        std_0_raw = rep(sqrt(my_var_raw),H),
        std_0_cen = rep(sqrt(my_var_cen),H))
}


#' Spatial Median
#'
#' Compute Spatial median of data
#'
#' @param tt Numeric for fparam (intraday) times.
#' @param x Matrix for data.
#' @param dtyp String 'n' or 's'.
#'
#' @returns Spatial median estimate.
#'
#' @noRd
#' @keywords internal
.spatial_median <- function(tt, x, dtyp = 's'){
  tt <- as.matrix(tt, ncol = length(tt), nrow = 1)
  eps <- 2^(-52)

  n <- nrow(x)
  m <- ncol(x)

  if (nrow(tt) > 1)
    tt <- t(tt)

  if (ncol(tt) != m)
    stop('Dimensions of T and X not compatible')

  if (!(dtyp != 'n' | dtyp != 's'))
    stop('Wrong input value for DTYP')

  # Inner-product matrix
  A <- 0.5 * (x[, 1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 1:(m-1)]) +
                x[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 2:m]))

  if( tolower(dtyp) == "n"){
    se2 <- .gsj(kronecker(matrix(1, 1, n), tt), t(x))
    A1 <- A - diag(se2) * (tt[m] - tt[1]);
    if (min(eigen(A1)>-1e-6))
      A <- A1
    else{
      message('Corrected inner-product matrix is not nonnegative definite \n Using uncorrected matrix instead ')
    }
  }

  ## Iterative minimization from sample mean
  w <- matrix(1, nrow = n, ncol = 1)/n # a very naive way
  norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # pay attention to the R syntax that it does not like to do arithmetic on a number to vector
  f <- sum(norms)
  err <- 1
  iter <- 0

  while( err > 1E-5 && iter < 50){
    iter <- iter + 1
    f0 <- f
    if (min(norms < eps)){
      i0 <- utils::find(norms < eps)
      w <- matrix(0, nrow = n, ncol = 1)
      w[i0] <- 1/length(i0)
    } else{
      w <- 1/norms
      w <- w/sum(w)
    }

    norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # same here
    f <- sum(norms)
    err <- abs(f/f0 - 1)
  }

  med <- t(w) %*% x

  list(med = med, w = w, norms = norms)
}


#' Temporary file for SACF
#'
#' @param x TODO
#' @param y TODO
#'
#' @returns TODO
#'
#' @noRd
#' @keywords internal
.gsj <- function(x, y){
  n <- nrow(x)
  m <- ncol(x)
  if (n == 1){
    x <-  t(x)
    y <- t(y)
    n <- nrow(x)
    m <- ncol(x)
  }

  if (n < 3)
    return(NULL)
  ind <- apply(x, 2, order)
  x <- apply(x, 2, sort)

  # Also sort the y accordingly using the indicator above
  for (j in 1:m) {
    y[ ,j] <- y[ind[ ,j], j]
  }

  a <- (x[3:n, ] - x[2:(n-1), ]) / (x[3:n, ] - x[1:(n-2), ])
  b <- (x[2:(n-1), ] - x[1:(n-2), ])/(x[3:n, ]-x[1:(n-2), ])
  c <- a^2 + b^2 + 1
  e <- a * y[1:(n-2), ] + b * y[3:n, ] - y[2:(n-1),]

  colMeans(e^2/c)
}


#' Temp Calculate Variance
#'
#' @param xbr TODO
#'
#' @returns TODO
#'
#' @noRd
#' @keywords internal
.calc_var <- function(xbr){
  Nt <- ncol(xbr)
  N <- nrow(xbr)
  res <- matrix(NA, Nt, Nt)
  for(s in 1:Nt){
    for(t in 1:Nt){
      val <- 0
      for(i in 1:N){
        val <- val + xbr[i, s] * xbr[i, t]
      }
      res[s, t] <- (val/N)^2
    }
  }
  sum(res)/(Nt-1)^2
}


#' Temp Calculate Norm
#'
#' @param obs TODO
#'
#' @returns TODO
#'
#' @noRd
#' @keywords internal
.calc_norm_matrix <- function(obs){
  Nt <- ncol(obs)
  N <- nrow(obs)
  sapply(1:N, function(i){
    sqrt(sum(obs[i, ]^2))/sqrt(Nt-1)
  })
}


#' Temp Calculate Variance Raw
#'
#' @param x_raw TODO
#'
#' @returns TODO
#'
#' @noRd
#' @keywords internal
.calc_var_raw <- function(x_raw){
  N <- nrow(x_raw)
  Nt <- ncol(x_raw)

  Spatial_med <- as.numeric(.spatial_median(tt = seq(0, 1, length.out = Nt), x_raw)$med)
  x_cen <- sapply(1:Nt, function(k){
    x_raw[, k] - Spatial_med[k]
  })

  my_norm <- .calc_norm_matrix(x_cen)
  xbr <- x_cen/my_norm

  res <- matrix(NA, Nt, Nt)
  for(s in 1:Nt){
    for(t in 1:Nt){
      val <- 0
      for(i in 1:N){
        val <- val + xbr[i, s] * xbr[i, t]
      }
      res[s, t] <- (val/N)^2
    }
  }
  sum(res)/(Nt-1)^2
}

