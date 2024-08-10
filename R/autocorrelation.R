#' Estimate the autocorrelation function of the series
#'
#' @description Obtain the empirical autocorrelation function for
#' lags \eqn{= 0,...,}\code{nlags} of the functional time
#' series. Given \eqn{Y_{1},...,Y_{T}} a functional time
#' series, the sample autocovariance functions
#' \eqn{\hat{C}_{h}(u,v)} are given by:
#' \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
#' where
#' \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
#' denotes the sample mean function. By normalizing these
#' functions using the normalizing factor
#' \eqn{\int\hat{C}_{0}(u,u)du}, the range of the
#' autocovariance functions becomes \eqn{(0,1)}; thus
#' defining the autocorrelation functions of the series
#'
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Discretization points of the curves, by default
#' \code{seq(from = 0, to = 1, length.out = 100)}.
#' @param nlags Number of lagged covariance operators
#' of the functional time series that will be used
#' to estimate the autocorrelation function.
#'
#' @return Return a list with the lagged autocorrelation
#' functions estimated from the data. Each function is given
#' by a \eqn{(m x m)} matrix, where \eqn{m} is the
#' number of points observed in each curve.
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' bbridge <- generate_brownian_bridge(N, v, sig)
#' nlags <- 1
#' lagged_autocor <- obtain_autocorrelation(Y = bbridge,
#'                                         nlags = nlags)
#' image(x = v, y = v, z = lagged_autocor$Lag0)
#'
#' \donttest{
#' # Example 2
#' require(fields)
#' N <- 500
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 4
#' lagged_autocov <- obtain_autocovariance(Y = bbridge,
#'                                         nlags = nlags)
#' lagged_autocor <- obtain_autocorrelation(Y = bbridge,
#'                                          v = v,
#'                                          nlags = nlags)
#'
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,2))
#' z_lims <- range(lagged_autocov$Lag1)
#' colors <- heat.colors(12)
#' image.plot(x = v,
#'            y = v,
#'            z = lagged_autocov$Lag1,
#'            legend.width = 2,
#'            zlim = z_lims,
#'            col = colors,
#'            xlab = "u",
#'            ylab = "v",
#'            main = "Autocovariance")
#' z_lims <- range(lagged_autocor$Lag1)
#' image.plot(x = v,
#'            y = v,
#'            z = lagged_autocor$Lag1,
#'            legend.width = 2,
#'            zlim = z_lims,
#'            col = colors,
#'            xlab = "u",
#'            ylab = "v",
#'            main = "Autocorrelation")
#' par(opar)
#' }
obtain_autocorrelation <- function(Y, nlags){

  fun.autocovariance <- .compute_autocovariance(Y,nlags)
  normalization.value <- dot_integrate_uneven(r = Y$intraobs, v = diag(fun.autocovariance$Lag0))
  fun.autocorrelation <- list()
  for(kk in 0:nlags){
    fun.autocorrelation[[paste("Lag",kk,sep = "")]] <-
      fun.autocovariance[[paste("Lag",kk,sep = "")]] / normalization.value
  }

  fun.autocorrelation
}



#' Generate a 3D plot of the autocovariance surface of a given FTS
#'
#' @description Obtain a 3D plot of the autocovariance surfaces of a
#' given functional time series. This visualization is
#' useful to detect any kind of dependency between
#' the discretization points of the series.
#'
#' @param fun.autocovariance A list obtained by
#' calling the function \code{obtain_autocovariance}.
#' @param lag An integer between 0 and \code{nlags}, indicating
#' the lagged autocovariance function to be plotted.
#' By default 0.
#' @param ... Further arguments passed to the  \code{persp}
#' function.
#'
#' @return XX
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' bbridge <- generate_brownian_bridge(N, v, sig)
#' nlags <- 1
#' lagged_autocov <- .compute_autocovariance(X = bbridge,nlags = nlags)
#' plot_autocovariance(lagged_autocov,1)
#'
#' \donttest{
#' # Example 2
#'
#' N <- 500
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 4
#' lagged_autocov <- obtain_autocovariance(Y = bbridge,nlags = nlags)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,5))
#' for(k in 0:nlags){
#'    plot_autocovariance(lagged_autocov,k)
#' }
#' par(opar)
#' }
plot_autocovariance <- function(fun.autocovariance, lag = 0, ...){

  # Color palette
  col.pal<-grDevices::colorRampPalette(c("blue", "red"))
  colors<-col.pal(100)

  # Select the lagged autocovariance surface
  z <- fun.autocovariance[[paste("Lag",lag,sep="")]]
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))     arguments$xlab <- "u"
  if(!"ylab" %in% names(arguments))     arguments$ylab <- "v"
  if(!"zlab" %in% names(arguments))     arguments$zlab <- ""
  if(!"main" %in% names(arguments))     arguments$main  <- paste("Lag",lag,sep=" ")
  if(!"expand" %in% names(arguments))   arguments$expand  <- 0.7
  if(!"ticktype" %in% names(arguments)) arguments$ticktype  <-'detailed'
  if(!"shade" %in% names(arguments))    arguments$shade  <- NA
  if(!"col" %in% names(arguments))      arguments$col <- colors[z.facet.range]
  if(!"theta" %in% names(arguments))    arguments$theta  <- 315
  if(!"phi" %in% names(arguments))      arguments$phi  <- 30
  arguments$z = z

  do.call(graphics::persp,arguments)
}

