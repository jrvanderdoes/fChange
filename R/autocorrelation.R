#' Estimate the autocorrelation function of the series
#'
#' Obtain the empirical autocorrelation function for lags
#'  \eqn{= 0,...,}\code{nlags} of the functional time series. Given
#'  \eqn{Y_{1},...,Y_{T}} a functional time series, the sample autocovariance
#'  functions \eqn{\hat{C}_{h}(u,v)} are given by:
#'  \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{X}_{T}(u))(Y_{i+h}(v) - \overline{X}_{T}(v))}
#'  where
#'  \eqn{ \overline{X}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
#'  denotes the sample mean function. By normalizing these functions using
#'  the normalizing factor \eqn{\int\hat{C}_{0}(u,u)du}, the range of the
#'  autocovariance functions becomes \eqn{(0,1)}; thus defining the
#'  autocorrelation functions of the series.
#'
#' @param X funts object
#' @param nlags Number of lagged covariance operators of the functional
#'  time series that will be used to estimate the autocorrelation function.
#'
#' @return Return a list with the lagged autocorrelation functions estimated
#'  from the data. Each function is given by a \eqn{(m x m)} matrix, where
#'  \eqn{m} is the number of points observed in each curve.
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
#' lagged_autocor <- obtain_autocorrelation(X = bbridge,
#'                                         nlags = nlags)
#' image(x = v, y = v, z = lagged_autocor$Lag0)
#'
#' \donttest{
#' # Example 2
#' require(fields)
#' N <- 500
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' bbridge <- generate_brownian_bridge(N, v, sig)
#' nlags <- 4
#' lagged_autocov <- obtain_autocovariance(X = bbridge,
#'                                         nlags = nlags)
#' lagged_autocor <- obtain_autocorrelation(X = bbridge,
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
obtain_autocorrelation <- function(X, nlags){
  X <- .check_data(X)

  fun.autocovariance <- .compute_autocovariance(X,nlags)
  normalization.value <-
    dot_integrate_uneven(r = X$intraobs, v = diag(fun.autocovariance$Lag0))

  fun.autocorrelation <- list()
  for(kk in 0:nlags){
    fun.autocorrelation[[paste("Lag",kk,sep = "")]] <-
      fun.autocovariance[[paste("Lag",kk,sep = "")]] / normalization.value
  }

  fun.autocorrelation
}


#' Generate a 3D plot of the autocovariance surface of a given FTS
#'
#' Obtain a 3D plot of the autocovariance surfaces of a given functional
#'  time series. This visualization is useful to detect any kind of dependency
#'  between the discretization points of the series.
#'
#' @param fun.autocovariance A list obtained by calling the function
#'  \code{obtain_autocovariance}.
#' @param lag An integer between 0 and \code{nlags}, indicating the lagged
#'  autocovariance function to be plotted. By default 0.
#' @param ... Further arguments passed to the  \code{persp} function.
#'
#' @return Graphic from \code{persp}
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
#' lagged_autocov <- obtain_autocovariance(X = bbridge,nlags = nlags)
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



#' Estimate the autocovariance function of the series
#'
#' Obtain the empirical autocovariance function for lags
#'   \eqn{= 0,...,}\code{nlags} of the functional time series. Given
#'   \eqn{Y_{1},...,Y_{T}} a functional time series, the sample autocovariance
#'   functions \eqn{\hat{C}_{h}(u,v)} are given by:
#'   \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
#'   where
#'   \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
#'   denotes the sample mean function.
#'
#' @param X funts object
#' @param nlags Number of lagged covariance operators of the functional
#'  time series that will be used to estimate the autocorrelation function.
#'
#' @return Return a list with the lagged autocovariance functions
#'  estimated from the data. Each function is given by a \eqn{(m x m)}
#'  matrix, where \eqn{m} is the number of points observed in each curve.
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 1
#' lagged_autocov <- obtain_autocovariance(Y = bbridge,
#'                                         nlags = nlags)
#' image(x = v, y = v, z = lagged_autocov$Lag0)
#'
#' \donttest{
#' # Example 2
#'
#' N <- 500
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 10
#' lagged_autocov <- obtain_autocovariance(Y = bbridge,
#'                                         nlags = nlags)
#' image(x = v, y = v, z = lagged_autocov$Lag0)
#' image(x = v, y = v, z = lagged_autocov$Lag10)
#'
#' # Example 3
#'
#' require(fields)
#' N <- 500
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 4
#' lagged_autocov <- obtain_autocovariance(Y = bbridge,
#'                                         nlags = nlags)
#' z_lims <- range(lagged_autocov$Lag0)
#' colors <- heat.colors(12)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,5))
#' par(oma=c( 0,0,0,6))
#' for(k in 0:nlags){
#'    image(x=v,
#'          y=v,
#'          z = lagged_autocov[[paste0("Lag",k)]],
#'          main = paste("Lag",k),
#'          col = colors,
#'          xlab = "u",
#'          ylab = "v")
#' }
#' par(oma=c( 0,0,0,2.5)) # reset margin to be much smaller.
#' image.plot( legend.only=TRUE, legend.width = 2,zlim=z_lims, col = colors)
#' par(opar)
#' }
obtain_autocovariance <- function(X, nlags){
  .compute_autocovariance(X, nlags)
}
