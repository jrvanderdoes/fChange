#' Kernel Functions
#'
#' There are an assortment of (vectorized) kernel functions located in the package.
#'
#' @param x A numeric value at which to evaluate kernel.
#'  It often indicates current lag divided by window.
#'
#' @name kernels
#'
#' @references Horvath, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @references L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional
#'  time series", Journal of Econometrics, 179(1), 66-82.
#'
#' @references Politis, D. N. (2003). Adaptive bandwidth choice. Journal of Nonparametric
#'  Statistics, 15(4-5), 517-533. \url{https://doi.org/10.1080/10485250310001604659}
#'
#' @references Politis, D. N. (2011). HIGHER-ORDER ACCURATE, POSITIVE SEMIDEFINITE
#'  ESTIMATION OF LARGE-SAMPLE COVARIANCE AND SPECTRAL DENSITY MATRICES.
#'  Econometric Theory, 27(4), 703-744. \url{http://www.jstor.org/stable/27975501}
#'
#' @return Values from a given lag in the kernel
NULL


#' *Truncated Kernel*: Kernel where \eqn{1, |x|\leq 1} and \eqn{0} otherwise. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' truncated_kernel(-20:20/15)
truncated_kernel <- function(x) {
  ifelse(is.nan(abs(x)), 1, ifelse(abs(x)<=1,1,0) )
}


#' *Bartlett Kernel*: Kernel where \eqn{max(0,1-|x|), h\neq 0}. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' bartlett_kernel(-20:20/15)
bartlett_kernel <- function(x) {
  pmax(0, ifelse(is.nan(1 - abs(x)),1,1 - abs(x)) )
}


#' *Parzen Kernel*: Kernel where \eqn{1 - 6 * x^2 + 6 * |x|^3, |x|<=0.5},
#'  \eqn{ 2 * (1 - |x|)^3, 0.5<|x|<1}, and \eqn{0, |x|>1}. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' parzen_kernel(-20:20/15)
parzen_kernel <- function(x) {
  ifelse(is.nan(x), 1,
    ifelse(abs(x) <= 1,
           ifelse(abs(x) <= 0.5,
                  1 - 6 * x^2 + 6 * abs(x)^3,
                  2 * (1 - abs(x))^3 ),
           0 )
  )
}


#' *Tukey-Hanning Kernel*: Kernel where \eqn{(1 + cos(\pi x) )/2, |x|<=1} and \eqn{0, |x|>1}. If
#'  \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' tukey_hanning_kernel(-20:20/15)
tukey_hanning_kernel <- function(x) {
  ifelse(is.nan(x), 1,
    ifelse(abs(x) <= 1,
           (1+cos(pi*x))/2,
           0 )
  )
}


#' *Quadratic Spectral Kernel*: Kernel where
#'  \eqn{\frac{25}{12\pi^2x^2} \left(\frac{sin(6\pi x/5)}{6\pi x/5} - cos(6\pi x/5) \right)}.
#'  If \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' quadratic_spectral_kernel(-20:20/15)
quadratic_spectral_kernel <- function(x) {
  res <- 25/(12 * pi^2 * x^2) * (sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5) )
  ifelse(is.nan(res), 1, res)
}


#' *Daniell Kernel*: Kernel where \eqn{sin(pi * x) / (pi * x)*(1 + cos(pi*x)), abs(x)<=1}.
#'  If \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' daniell_kernel(-20:20/15)
daniell_kernel <- function(x) {
  res <- ifelse(is.nan(x), 1,
                ifelse(abs(x) <= 1, sin(pi * x) / (pi * x) * (1 + cos(pi*x)), 0 )
  )
  ifelse( is.nan(res), 1, res )
}


#' *Flat-Top Kernel*: Kernel where \eqn{min(1, max(1.1-|x|,0)),|x|\leq 1}. If
#'  \eqn{x=0/0} then the value \eqn{1} is given. Note, values rather than 1.1
#'  have also been used.
#'
#' @rdname kernels
#'
#' @export
#'
#' @examples
#' flat_top_kernel(-20:20/15)
flat_top_kernel <- function(x) {
  # TODO:: Setup functions to allow FT with arbitrary width
  ifelse(is.nan(x), 1, pmin(1, pmax(1.1 - abs(x), 0)) )
}


#' Adaptive_bandwidth
#'
#' Computes the "optimal" bandwidth using a bandwidth selection method based on the
#'   spectral density operator which adapts to the functional data.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param kernel Kernel function selection to use, 'bartlett' or 'parzen'.
#'
#' @return a scalar value of the "optimal" data-adapted bandwidth.
#' @export
adaptive_bandwidth <- function(X, kernel='bartlett') {
  X <- dfts(X)
  f_data <- X$data

  ## TODO:: Setup to not be hidden
  J <- NROW(f_data)
  N <- NCOL(f_data)

  if (tolower(kernel) == 'bartlett') {
    kernel_fun <- bartlett_kernel
    order <- 1
    xi <- 1
    kern_int <- 2 / 3
  } else if (tolower(kernel) == 'parzen') {
    kernel_fun <- parzen_kernel
    order <- 2
    xi <- 6
    kern_int <- 151 / 280
    # } else if (kernel == 'Daniell') {
    # kernel_fun <- daniell_kernel
    # order <- 2
    # xi <- (pi^2) / 6
    # kern_int <- 1
  } else {
    stop('Please see the documentation for supported kernels.')
  }

  data_inner_prod <- crossprod(f_data) / (N * J)
  C_hat_HS <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS[j+1] <- sum( data_inner_prod[(j+1):N,(j+1):N] *
                          data_inner_prod[1:(N-j),1:(N-j)] )
  }

  initial_band_q <- 4 * N^(1 / (2*order +1))
  k_n_j_q <- kernel_fun(1:(N-1) / initial_band_q)
  initial_band_0 <- initial_band_q / 4
  k_n_j_0 <- kernel_fun(1:(N-1) / initial_band_0)

  Q_hat_sq <- 2 * sum(k_n_j_0^2 * C_hat_HS[-1])
  Q_hat_sq <- Q_hat_sq + sum(C_hat_HS[0])

  Term1 <- 2 * sum(k_n_j_q^2 * ((1:(N-1))^(2*order)) * C_hat_HS[-1])
  Term2 <- Q_hat_sq
  C_hat_TR <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_TR[j+1] <- sum(data_inner_prod[1:(N-j), (j+1):N])
  }
  trace <- 2 * sum(k_n_j_0^2 * C_hat_TR[-1])
  Term3 <- trace + sum(C_hat_TR[1])
  band_constant <- (2 * order * xi^2 * Term1 / (kern_int * (Term2 + Term3)))^(1/(2*order + 1))

  band_constant * N^(1 / (2*order + 1))
}

