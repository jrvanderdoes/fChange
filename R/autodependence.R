#' Estimate the autocovariance function of the series
#'
#' Obtain the empirical autocovariance function for the given \code{lags} of a
#'  functional time series, \code{X}. Given a functional time series, the sample
#'  autocovariance functions \eqn{\hat{C}_{h}(u,v)} are given by:
#'    \deqn{\hat{C}_{h}(u,v) =  \frac{1}{N} \sum_{i=1}^{N-|h|}(Y_{i}(u) -
#'      \overline{X}_{N}(u))(Y_{i+|h|}(v) - \overline{X}_{N}(v))}
#'  where \eqn{ \overline{X}_{N}(u) = \frac{1}{N} \sum_{i = 1}^{N} X_{i}(t)}
#'  denotes the sample mean function and \eqn{h} is the lag parameter.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param lags Numeric(s) for the lags to estimate the lagged operator.
#' @param center Boolean if the data should be centered. Default is true.
#'
#' @return Return a list or data.frame with the lagged autocovariance function(s)
#'  estimated from the data. Each function is given by a \eqn{(r } x \eqn{ r)}
#'  matrix, where \eqn{r} is the number of points observed in each curve.
#' @export
#'
#' @seealso [autocorrelation()], [var()]
#'
#' @examples
#' v <- seq(0,1,length.out=20)
#' lagged_autocov <- autocovariance(
#'   X = generate_brownian_bridge(100,v=v),
#'   lags = 1)
autocovariance <- function(X, lags=0:1, center=TRUE){
  X <- dfts(X)
  if(center) X <- center.dfts(X)
  res <- length(X$fparam)
  obs <- length(X$labels)

  autocovs <- list()
  for(nlag in lags){
    # autocovs[[paste0("Lag", nlag)]] <-
    #   matrix(0, nrow = res, ncol = res)
    idxs <- (1+nlag):obs

    autocovs[[paste0("Lag", nlag)]] <-
      X$data[,idxs - nlag,drop=FALSE] %*% t(X$data[,idxs,drop=FALSE]) /
      obs#(obs-1) #TODO:: n or n-1
  }

  if(length(lags)==1) return(autocovs[[1]])

  autocovs
}


#' Estimate the autocorrelation function of the series
#'
#' Obtain the empirical autocorrelation function for the given \code{lags} of a
#'  functional time series, \code{X}. Given a functional time series, the sample
#'  autocovariance functions \eqn{\hat{C}_{h}(u,v)} are given by:
#'    \deqn{\hat{C}_{h}(u,v) =  \frac{1}{N} \sum_{i=1}^{N-|h|}(X_{i}(u) -
#'      \overline{X}_{N}(u))(X_{i+|h|}(v) - \overline{X}_{N}(v))}
#'  where \eqn{ \overline{X}_{N}(u) = \frac{1}{N} \sum_{i = 1}^{N} X_{i}(t)}
#'  denotes the sample mean function and \eqn{h} is the lag parameter. The
#'  autocorrelation functions are defined over the range \eqn{(0,1)} by
#'  normalizing these functions using the factor \eqn{\int\hat{C}_{0}(u,u)du}.
#'
#' @inheritParams autocovariance
#'
#' @return Return a list or data.frame with the lagged autocorrelation function(s)
#'  estimated from the data. Each function is given by a \eqn{(r } x \eqn{ r)}
#'  matrix, where \eqn{r} is the number of points observed in each curve.
#' @export
#'
#' @seealso [autocovariance()]
#'
#' @examples
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' bbridge <- generate_brownian_bridge(N = N, v = v)
#' lagged_autocor <- autocorrelation(X = bbridge, lags = 0:1)
autocorrelation <- function(X, lags){
  X <- dfts(X)

  lags_use <- unique(c(0,lags))

  autocov <- autocovariance(X, lags_use)

  if(length(lags_use)==1){
    normalization.value <-
      dot_integrate(r = X$fparam, v = diag(autocov))

    autocor <- autocov / normalization.value

  } else if(length(lags)==1){
    normalization.value <-
      dot_integrate(r = X$fparam, v = diag(autocov$Lag0))

    autocor <- autocov[[paste("Lag",lags,sep = "")]] / normalization.value

  }else {
    normalization.value <-
      dot_integrate(r = X$fparam, v = diag(autocov$Lag0))

    autocor <- list()
    for(kk in lags){
      autocor[[paste("Lag",kk,sep = "")]] <-
        autocov[[paste("Lag",kk,sep = "")]] / normalization.value
    }
  }


  autocor
}
