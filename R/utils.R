#' Approximate Integral
#'
#' This (internal) function computes the estimates the integral for
#'  evenly spaced observations.
#'
#' @param y Numeric vector of data
#' @param type (Optional) String indicating type of integration. Currently only
#'  the only option is "Rectangle" for rectangular integration done using
#'  Reimann Sum. Default is "Rectangle".
#'
#' @return Numeric indicating the estimated integral of the curve
#'
#' @noRd
#'
#' @examples
#' .approx_int(rep(1,10))
#' .approx_int(seq(0,1,length.out=20))
#' .approx_int(seq(0,1,length.out=20)^2)
.approx_int <- function(y, type="Rectangle"){
  value <- 0
  # Trap Approx 1
  # pracma::trapz(seq(0,1,length.out=length(y)+1), c(0,y))
  # Trap Approx 2
  # pracma::trapz(seq(0,1,length.out=length(y)),y)
  # Trap Approx 3
  # n <- length(y)
  # pracma::trapz(seq(1/n,1,length.out=n),y)
  # Trap Approx 4
  # n <- length(y)
  # pracma::trapz(seq(0,1-1/n,length.out=n), y)
  if(type=="Rectangle"){ # Rect Approx
    value <- sum(y) / length(y)
    # 1/(length(y)+1)*sum(c(0,y))
  } else if(type=="Trapezoidal"){ # Trap Approx
    value <- sum(y[-1] + y[length(y)]) / ( 2 * length(y) )
    stop('Error: Trapezoidal rule not yet setup.')
    # pracma::trapz(seq(0,1,length.out=length(y)),y)
    # pracma::trapz(seq(0,1,length.out=length(y)),y)
  }

  value
}
