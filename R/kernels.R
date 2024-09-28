#' Truncated Kernel
#'
#' Kernel where \eqn{1, |x|\leq 1} and \eqn{0} otherwise. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @param x A numeric value at which to evaluate kernel.
#'  It often indicates current lag divided by window.
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @examples
#' truncated_kernel(-20:20/15)
truncated_kernel <- function(x) {
  ifelse(is.nan(abs(x)), 1, ifelse(abs(x)<=1,1,0) )
}


#' Bartlett Kernel
#'
#' Kernel where \eqn{max(0,1-|x|), h\neq 0}. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @param x A numeric value at which to evaluate kernel.
#'  It often indicates current lag divided by window.
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @examples
#' bartlett_kernel(-20:20/15)
bartlett_kernel <- function(x) {
  pmax(0, ifelse(is.nan(1 - abs(x)),1,1 - abs(x)) )
}


#' Parzen Kernel
#'
#' Kernel where \eqn{1 - 6 * x^2 + 6 * |x|^3, |x|<=0.5},
#'  \eqn{ 2 * (1 - |x|)^3, 0.5<|x|<1}, and \eqn{0, |x|>1}. If \eqn{x=0/0} then
#'  the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
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


#' Tukey-Hanning Kernel
#'
#' Kernel where \eqn{(1 + cos(\pi x) )/2, |x|<=1} and \eqn{0, |x|>1}. If
#'  \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
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


#' Quadratic Spectral Kernel
#'
#' Kernel where \eqn{\frac{25}{12\pi^2x^2} \left(\frac{sin(6\pi x/5)}{6\pi x/5} - cos(6\pi x/5) \right)}.
#'  If \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @examples
#' quadratic_spectral_kernel(-20:20/15)
quadratic_spectral_kernel <- function(x) {
  res <- 25/(12 * pi^2 * x^2) * (sin(6*pi*x/5)/(6*pi*x/5) - cos(6*pi*x/5) )
  ifelse(is.nan(res), 1, res)
}


#' Daniell Kernel
#'
#' Kernel where \eqn{sin(pi * x) / (pi * x)*(1 + cos(pi*x)), abs(x)<=1}.
#'  If \eqn{x=0/0} then the value \eqn{1} is given.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @examples
#' daniell_kernel(-20:20/15)
daniell_kernel <- function(x) {
  res <- ifelse(is.nan(x), 1,
                ifelse(abs(x) <= 1, sin(pi * x) / (pi * x) * (1 + cos(pi*x)), 0 )
  )
  ifelse( is.nan(res), 1, res )
}


#' Flat-Top Kernel
#'
#' Kernel where \eqn{min(1, max(1.1-|x|,0)),|x|\leq 1}. If \eqn{x=0/0} then
#'  the value \eqn{1} is given. Note, values rather than 1.1 have also been
#'  used.
#'
#' Note, this function is vectorized.
#'
#' @inheritParams bartlett_kernel
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#'
#' @references L. Horvath, P. Kokoszka, G. Rice (2014) "Testing stationarity of functional
#'  time series", Journal of Econometrics, 179(1), 66-82.
#' @references Politis, D. N. (2003). Adaptive bandwidth choice. Journal of Nonparametric
#'  Statistics, 15(4–5), 517–533. \url{https://doi.org/10.1080/10485250310001604659}
#' @references Politis, D. N. (2011). HIGHER-ORDER ACCURATE, POSITIVE SEMIDEFINITE
#'  ESTIMATION OF LARGE-SAMPLE COVARIANCE AND SPECTRAL DENSITY MATRICES.
#'  Econometric Theory, 27(4), 703–744. \url{http://www.jstor.org/stable/27975501}
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (1st ed. 2024.). Springer Nature Switzerland.
#'  \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @examples
#' flat_top_kernel(-20:20/15)
flat_top_kernel <- function(x) {
  # TODO:: Setup functions to allow FT with arbitrary width
  ifelse(is.nan(x), 1, pmin(1, pmax(1.1 - abs(x), 0)) )
}
