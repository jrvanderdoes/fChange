#' Kernel Functions
#'
#' There are an assortment of (vectorized) kernel functions located in the package.
#'
#' @param x Numeric value(s) at which to evaluate kernel.
#'  It often indicates current lag divided by window.
#'
#' @name kernels
#'
#' @references Horvath, L., & Rice, G. (2024). Change point analysis for time
#'  series (1st ed. 2024.). Springer Nature Switzerland.
#'
#' @references L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional
#'  time series", Journal of Econometrics, 179(1), 66-82.
#'
#' @references Politis, D. N. (2003). Adaptive bandwidth choice. Journal of Nonparametric
#'  Statistics, 15(4-5), 517-533.
#'
#' @references Politis, D. N. (2011). Higher-order accurate, positive semidefinite
#'  estimation of large-sample covariance and spectral density matrices.
#'  Econometric Theory, 27(4), 703-744.
#'
#' @return Values from given lag(s) in the kernel.
NULL


#' @description *Truncated Kernel*: Kernel where \eqn{1, |x|\leq 1} and \eqn{0}
#'  otherwise. If \eqn{x=0/0} then the value \eqn{1} is given.
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


#' @description *Bartlett Kernel*: Kernel where \eqn{max(0,1-|x|), h\neq 0}. If
#'  \eqn{x=0/0} then the value \eqn{1} is given.
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


#' @description *Parzen Kernel*: Kernel where \eqn{1 - 6 * x^2 + 6 * |x|^3, |x|<=0.5},
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


#' @description *Tukey-Hanning Kernel*: Kernel where \eqn{(1 + cos(\pi x) )/2, |x|<=1}
#'  and \eqn{0, |x|>1}. If \eqn{x=0/0} then the value \eqn{1} is given.
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


#' @description *Quadratic Spectral Kernel*: Kernel where
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


#' @description *Daniell Kernel*: Kernel where \eqn{sin(pi * x) / (pi * x)*(1 + cos(pi*x)), abs(x)<=1}.
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


#' @description *Flat-Top Kernel*: Kernel where \eqn{min(1, max(1.1-|x|,0)),|x|\leq 1}.
#'  If \eqn{x=0/0} then the value \eqn{1} is given.
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
#' Computes the data-driven bandwidth using a method based on the
#'   spectral density operator which adapts to the functional data.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param kernel Kernel function. No additional parameters are needed for
#'  \code{bartlett_kernel()}, \code{parzen_kernel()}, \code{tukey_hanning_kernel()},
#'  and \code{quadratic_spectral_kernel()}.
#' @param name,order,weighting Additional parameters if non-standard kernels, e.g.
#'  those not in fChange, are used. See references for the definitions. Name is
#'  extracted from the kernel name to select \code{order}/\code{weighting} when
#'  not given, if the function aligns with the recommended functions, see
#'  \code{kernel} parameter.
#'
#' @references Rice, G., & Shang, H. L. (2017). A Plug‐in Bandwidth Selection
#'  Procedure for Long‐Run Covariance Estimation with Stationary Functional Time
#'  Series. Journal of Time Series Analysis, 38(4), 591–609.
#'
#' @seealso [bartlett_kernel()], [truncated_kernel()], [parzen_kernel()],
#'  [tukey_hanning_kernel()], [quadratic_spectral_kernel()], [daniell_kernel()],
#'  [flat_top_kernel()]
#'
#' @return Scalar value of the data-adapted bandwidth.
#' @export
#'
#' @examples
#' adaptive_bandwidth(generate_brownian_motion(100))
#' adaptive_bandwidth(electricity, parzen_kernel)
adaptive_bandwidth <- function(X, kernel=bartlett_kernel,
                               name=NULL, order=NULL, weighting=NULL) {
  X <- dfts(X)

  r <- nrow(X$data)
  N <- ncol(X$data)

  # Get kernel name if not provided for details

  if(is.null(order) || is.null(weighting)){
    if(is.null(name)) name <- deparse(substitute(kernel))
    tmp <- .kernel_details(kernel, name)

    if(is.null(order)) order <- tmp$order
    if(is.null(weighting)) weighting <- tmp$weighting
  }

  vals_int <- seq(-100,100,0.0001)
  kern_int <- dot_integrate(kernel(vals_int)^2, vals_int)

  # Code: Paper
  # N: T
  # order: q
  # kernel: W
  data_inner_prod <- crossprod(X$data) / (N * r)
  C_hat_HS <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_HS[j+1] <- sum( data_inner_prod[(j+1):N,(j+1):N] *
                            data_inner_prod[1:(N-j),1:(N-j)] )
  }

  # Pilot info C^(p), (6)
  initial_band_q <- 4 * N^(1 / (2*order +1))
  k_n_j_q <- kernel(1:(N-1) / initial_band_q)
  initial_band_0 <- initial_band_q / 4
  k_n_j_0 <- kernel(1:(N-1) / initial_band_0)

  Q_hat_sq <- 2 * sum(k_n_j_0^2 * C_hat_HS[-1])
  Q_hat_sq <- Q_hat_sq + sum(C_hat_HS[0])

  Term1 <- 2 * sum(k_n_j_q^2 * ((1:(N-1))^(2*order)) * C_hat_HS[-1]) # C^(q)
  C_hat_TR <- numeric(0)
  for (j in 0:(N-1)) {
    C_hat_TR[j+1] <- sum(data_inner_prod[1:(N-j), (j+1):N])
  }
  trace <- 2 * sum(k_n_j_0^2 * C_hat_TR[-1])
  Term3 <- trace + sum(C_hat_TR[1])
  band_constant <- (2 * order * weighting^2 * Term1 /
                      ( (Q_hat_sq + Term3) * kern_int) )^(1/(2*order + 1)) # After (5)

  band_constant * N^(1 / (2*order + 1)) # (5)
}


.kernel_details <- function(kernel, name){
  if( grepl('bartlett', tolower(name)) ) {
    order <- 1
    weighting <- 1
  } else if( grepl('parzen',tolower(name)) ){
    order <- 2
    weighting <- 6
  } else if( grepl('tukey',tolower(name)) && grepl('hanning',tolower(name)) ){
    order <- 2
    weighting <- (pi^2) / 6
  } else if( grepl('quadratic',tolower(name)) ){
    order <- 2
    weighting <- (18*pi^2) / 125
  } else{
    ## TODO:: flat_top, truncated_kernel, daniell_kernel
    stop('kernel variable not named with "barlett", "parzen", "tukey_hanning",
         "quadratic". Rename function or define name, order, and/or weighting.')

  }

  list('order'=order, 'weighting'=weighting)
}
