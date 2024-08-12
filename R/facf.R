#' Obtain the autocorrelation function for a given functional time series.
#'
#' Estimate the lagged autocorrelation function for a given
#'  functional time series and its distribution under the
#'  hypothesis of strong functional white noise. This graphic tool
#'  can be used to identify seasonal patterns in the functional
#'  data as well as auto-regressive or moving average terms.
#'  Strong white noise is included and weak white noise can also be included to
#'  test the presence of serial correlation in the data.
#'
#' @param X funts object
#' @param lag.max Number of lagged covariance operators
#'  of the functional time series that will be used
#'  to estimate the autocorrelation function. Default of NULL results in 20 lags.
#' @param alpha A value between 0 and 1 that indicates significant level for
#'  the confidence interval for the i.i.d. bounds of the autocorrelation
#'  function. By default \code{alpha = 0.95}.
#' @param method Character specifying the method to be used when estimating the
#'  distribution under the hypothesis of functional white noise.
#'  Accepted values are:
#'  \itemize{
#'    \item "Welch": Welch approximation.
#'    \item "MC": Monte-Carlo estimation.
#'    \item "Imhof": Estimation using Imhof's method.
#'  }
#'  By default, \code{method = "Welch"}.
#' @param WWN Logical. If \code{TRUE}, WWN bounds are also computed
#' @param figure Logical. If \code{TRUE}, plots the estimated autocorrelation
#'  function with the specified bounds.
#' @param ... Further arguments passed to the \code{.plot_FACF} function.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{SWN_bound}: The upper prediction
#'     bound for the i.i.d. distribution under strong white noise assumption.
#'     \item \code{WWN_bound}: The upper prediction
#'     bound for the i.i.d. distribution under weak white noise assumption.
#'     \item \code{acfs}: Autocorrelation values for
#'     each lag of the functional time series.
#' }
#' @export
#'
#' @references
#' Mestre G., Portela J., Rice G., Muñoz San Roque A., Alonso E. (2021).
#'  \emph{Functional time series model identification and diagnosis by
#'  means of auto- and partial autocorrelation analysis.}
#'  Computational Statistics & Data Analysis, 155, 107108.
#'  \url{https://doi.org/10.1016/j.csda.2020.107108}
#'
#' Mestre, G., Portela, J., Muñoz-San Roque, A., Alonso, E. (2020).
#'  \emph{Forecasting hourly supply curves in the Italian Day-Ahead
#'  electricity market with a double-seasonal SARMAHX model.}
#'  International Journal of Electrical Power & Energy Systems,
#'  121, 106083. \url{https://doi.org/10.1016/j.ijepes.2020.106083}
#'
#' Kokoszka, P., Rice, G., Shang, H.L. (2017).
#'  \emph{Inference for the autocovariance of a functional
#'  time series under conditional heteroscedasticity}
#'  Journal of Multivariate Analysis,
#'  162, 32--50. \url{https://doi.org/10.1016/j.jmva.2017.08.004}
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 5)
#' sd <- 2
#' X <- generate_brownian_bridge(N, v, sd)
#' .compute_FACF(X,20)
#'
#' \donttest{
#' # Example 2
#'
#' .compute_FACF(X=funts(electricity),lag.max = 30,alpha = 0.90)
#' }
.compute_FACF <- function(X, lag.max = NULL, alpha=0.05,
                          method = c('Welch','MC','Imhof'),
                          WWN = TRUE, figure = TRUE, ...){
  X <- .check_data(X)

  if(is.null(lag.max))
    lag.max <- 20 # 10 * log(ncol(X)/, base=10)

  # Autocov surfaces
  autocovs <- .compute_autocovariance(X, lag.max)

  # L2 norm autocov surfaces
  l2norms <- .obtain_suface_L2_norm(X$intraobs, autocovs)
  l2norms <- l2norms[-1] # Drop Lag 0

  # Obtain autocorrelation estimates
  normalization.value <-
    dot_integrate_uneven(r = X$intraobs, v = diag(autocovs$Lag0))
  rho <- sqrt(l2norms) / normalization.value

  # Estimate distribution of SWN (iid) bound
  method <- match.arg(method, c('Welch','MC','Imhof'))
  if(method=='Welch'){
    # Obtain SWN bound for specified confidence value
    SWN_bound <- rep(NA,length(alpha))
    for(ii in 1:length(alpha)){
      SWN_bound[ii] <-
        sqrt(Q_WS_quantile_iid(X$data, alpha=alpha)$quantile) /
        ( sqrt(NCOL(X$data)) * normalization.value )
    }
  } else{
    if(method=='MC'){
      iid.distribution <- .estimate_iid_distr_MC(X, autocovs, l2norms)
    } else if(method=='Imhof'){
      iid.distribution <- .estimate_iid_distr_Imhof(X, autocovs, l2norms)
    }

    # Obtain SWN bound for specified confidence value
    SWN_bound <- rep(NA,length(alpha))
    for(ii in 1:length(alpha)){
      idx <- min(which(iid.distribution$ef >= 1-alpha[ii]))
      SWN_bound[ii] <-
        sqrt( iid.distribution$ex[ idx ] ) /
        normalization.value
    }
  }

  # Obtain WWN bound for specified confidence value
  WWN_bound <- array(NA, lag.max)
  lags <- 1:lag.max

  if(WWN){
    quantile <- Q_WS_quantile(X$data, lags, alpha=alpha, M=NULL, low_disc=FALSE)$quantile
    WWN_bound <- sqrt(quantile) /
      ( sqrt(ncol(X$data)) * normalization.value )
  }

  # Plot
  if(figure){
    .plot_FACF(rho,SWN=SWN_bound, WWN=WWN_bound, ...)
  }

  invisible( list('SWN_bound'=SWN_bound, 'WWN_bound'=WWN_bound, 'acfs'=rho) )
}


#' Estimate Autocovariance Function of funts
#'
#' Obtain the empirical autocovariance function for lags
#'  \eqn{= 0,...,}\code{lag.max} of the functional time series. Given
#'  \eqn{X_{1},...,X_{T}} a functional time series, the sample
#'  autocovariance functions \eqn{\hat{C}_{h}(u,v)} are given by:
#'  \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
#'  where
#'  \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
#'  denotes the sample mean function.
#'
#' @inheritParams .compute_FACF
#'
#' @return Return a list with the lagged autocovariance
#' functions estimated from the data. Each function is given
#' by a \eqn{(m x m)} matrix, where \eqn{m} is the
#' resolution observed in each curve.
#' @export
#'
#' @examples
#' .compute_autocovariance(funts(electricity), 10)
.compute_autocovariance <- function(X, lag.max){
  X <- .check_data(X)
  lag.max <- ifelse(is.null(lag.max),20,lag.max)

  X_demean <- X$data - rowMeans(X$data) # TODO:: Check row means
  res <- nrow(X$data)
  obs <- ncol(X$data)

  autocovs <- list()
  for(nlag in 0:lag.max){
    autocovs[[paste0("Lag", nlag)]] <-
      matrix(0, nrow = res, ncol = res)
    for( ind_curve in (1+nlag):obs ){
      autocovs[[paste0("Lag", nlag)]] <-
        autocovs[[paste0("Lag", nlag)]] +
        X_demean[,ind_curve - nlag] %*% t(X_demean[,ind_curve])
    }
    autocovs[[paste0("Lag", nlag)]] <-
      autocovs[[paste0("Lag", nlag)]] / obs #TODO:: n or n-1
  }

  autocovs
}


#' Obtain L2 norm of the autocovariance functions
#'
#' Returns the L2 norm of the lagged autocovariance functions
#'  \eqn{\hat{C}_{h}}. The L2 norm of these functions is defined as
#'  \deqn{\sqrt(\int \int \hat{C}^{2}_{h}(u,v)du dv)}.
#'
#' @param intraobs Discretization points of the curves.
#' @param autocovs An \eqn{(m x m)} matrix with the discretized
#'  values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained by
#'  calling the function \code{obtain_autocovariance}. The value
#'  \eqn{m} indicates the number of points observed in each curve.
#'
#' @return A vector containing the L2 norm of the
#' lagged autocovariance functions \code{autocovs}.
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sd <- 2
#' X <- generate_brownian_bridge(N, v, sd)
#' lag.max <- 1
#' autocovSurface <- obtain_autocovariance(X=X,nlags = lag.max)
#' norms <- .obtain_suface_L2_norm(intraobs = v,autocovs = autocovSurface)
#' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 1)
#' title(sub = paste0("Lag ",1," - L2 Norm: ",norms[2]))
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' X <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 2
#' autocovSurface <- obtain_autocovariance(X=X,nlags = lag.max)
#' norms <- .obtain_suface_L2_norm(intraobs = v,autocovs = autocovSurface)
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,3))
#' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 0)
#' title(sub = paste0("Lag ",0," - L2 Norm: ",norms[1]))
#' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 1)
#' title(sub = paste0("Lag ",1," - L2 Norm: ",norms[2]))
#' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 2)
#' title(sub = paste0("Lag ",2," - L2 Norm: ",norms[3]))
#' par(opar)
#' }
.obtain_suface_L2_norm <- function(intraobs, autocovs){
  lag.max <- length(autocovs)
  matindex <- rep(NA, lag.max)
  res <- nrow(autocovs$Lag0)

  # L2 norm of squared autocov surfaces
  for(rr in 1:lag.max){
    autocov2 <- autocovs[[paste("Lag",rr-1,sep="")]]^2

    norm.vec <- matrix(NA, nrow = res, ncol = 1)
    for(ii in 1:res){
      norm.vec[ii] <- dot_integrate_uneven(r = intraobs, v = autocov2[,ii])
    }

    norm.aux <- dot_integrate_uneven(r = intraobs, v = norm.vec)
    matindex[rr] <- norm.aux
  }

  matindex
}


#' Estimate distribution of the fACF under the iid hypothesis using MC method
#'
#' Estimate the distribution of the autocorrelation function under the
#'  hypothesis of strong functional white noise. This function uses a
#'  Monte Carlo method to estimate the distribution.
#'
#' @inheritParams .compute_FACF
#' @param autocovSurface An \eqn{(m x m)} matrix with the discretized
#' values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained
#' by calling the function \code{obtain_autocovariance}.
#' The value \eqn{m} indicates the number of points observed
#' in each curve.
#' @param matindex A vector containing the L2 norm of
#' the autocovariance function. It can be obtained by calling
#' function \code{obtain_suface_L2_norm}.
#' @param nsims Positive integer indicating the number of
#' MC simulations that will be used to estimate the distribution
#' of the statistic. Increasing the number of simulations will
#' improve the estimation, but it will increase the computational
#' time.
#' By default, \code{nsims = 10000}.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{ex}: Knots where the distribution has been estimated.
#'     \item \code{ef}: Discretized values of the estimated distribution.
#'     \item \code{Reig}: Raw values for iid statistic of each MC simulation.
#' }
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 1
#' autocovSurface <- obtain_autocovariance(Y,lag.max)
#' matindex <- .obtain_suface_L2_norm(v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' MC_dist <- .estimate_iid_distr_MC(Y,autocovSurface,matindex)
#' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
#' grid()
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' Y <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 20
#' autocovSurface <- obtain_autocovariance(Y,lag.max)
#' matindex <- .obtain_suface_L2_norm(v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' MC_dist <- .estimate_iid_distr_MC(Y,autocovSurface,matindex)
#' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
#' grid()
#' }
.estimate_iid_distr_MC <-
  function(X, autocovSurface, matindex, nsims= 10000){
    X <- .check_data(X)

    # TODO:: Update Means
    # # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
    # # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
    # means <- matrix(rep(rowMeans(X$data),ncol(X$data)), nrow=nrow(X$data))
    # l <- .obtain_autocov_eigenvalues(X$data - means,X$intraobs)
    l <- .obtain_autocov_eigenvalues(center(X))

    neig <- length(l)
    Reig <- rep(0,nsims)

    for(jj in 1:neig){
      for(kk in 1:neig){
        Reig <- Reig +
          stats::rchisq(n = nsims,df = 1)*l[jj]*l[kk]
      }
    }
    #Reig=Reig/nrow(Y)
    Reig <- Reig / ncol(X$data)

    ecdf.aux <- stats::ecdf(Reig)
    ex <- stats::knots(ecdf.aux)
    ef <- ecdf.aux(ex)
    # if(figure){
    #   # Check if any additional plotting parameters are present
    #   arguments <- list(...)
    #   if("main" %in% names(arguments)){
    #     graphics::plot(ex,ef,type = "l",...)
    #   }else{
    #     graphics::plot(ex,ef,type = "l",main = "ecdf obtained by MC simulation",...)
    #   }
    #   graphics::grid()
    # }

    list('ex'=ex,'ef'=ef, 'Reig'=Reig)
  }


#' Estimate eigenvalues of the autocovariance function
#'
#' @description Estimate the eigenvalues of the sample autocovariance
#'  function \eqn{\hat{C}_{0}}. This functions returns the eigenvalues which
#'  are greater than the value \code{epsilon}.
#'
#' @inheritParams .compute_FACF
#' @param epsilon Value used to determine how many eigenvalues will be returned.
#'   The eigenvalues \eqn{\lambda_{j}/\lambda_{1} > \code{epsilon}} will be
#'   returned. By default \code{epsilon = 0.0001}.
#'
#' @return A vector containing the \eqn{k} eigenvalues
#' greater than \code{epsilon}.
#' @export
#'
#' @examples
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- generate_brownian_bridge(N, v, sig)
#' lambda <- .obtain_autocov_eigenvalues(X = Y)
.obtain_autocov_eigenvalues <- function(X, epsilon = 0.0001){
  X <- .check_data(X)

  nobs <- ncol(X$data) # nt
  res <- nrow(X$data) # nv

  # w(i,j) = integral of the product of the curves i and j
  W <- matrix(0,nrow = nobs,ncol = nobs)

  for(ii in 1:nobs){
    mat.aux <- matrix(rep(X$data[,ii],each = nobs),
                      nrow = nobs, ncol = res)*t(X$data)
    for(jj in 1:nobs){
      W[ii,jj] <- dot_integrate_uneven(r = X$intraobs, v = mat.aux[jj,])
    }
  }

  # Obtain eigenvalues and eigenfunctions of W/n
  eigenlist <- eigen(W/nobs)

  # Obtain all first $m$ eigenvalues such that
  #   \lambda_{m}/\lambda_{1} > \varepsilon
  lambd1 <- eigenlist$values[1]
  frac <- eigenlist$values/lambd1
  k <- which(frac > epsilon)

  eigenlist$values[k]
}


#' Estimate distribution of the fACF under the iid. hypothesis using Imhof's method
#'
#' Estimate the distribution of the autocorrelation function
#' under the hypothesis of strong functional white noise. This
#' function uses Imhof's method to estimate the distribution.
#'
#' @inheritParams .obtain_suface_L2_norm
#' @param l2norms A vector containing the L2 norm of
#' the autocovariance function. It can be obtained by calling
#' function \code{obtain_suface_L2_norm}.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{ex}: Knots where the distribution has been estimated.
#'     \item \code{ef}: Discretized values of the estimated distribution.
#' }
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 1
#' autocovSurface <- obtain_autocovariance(Y,lag.max)
#' matindex <- .obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' Imhof_dist <- .estimate_iid_distr_Imhof(Y,autocovSurface,matindex)
#' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
#' grid()
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' Y <- generate_brownian_bridge(N, v, sig)
#' autocovSurface <- obtain_autocovariance(Y,lag.max)
#' matindex <- .obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' Imhof_dist <- .estimate_iid_distr_Imhof(Y,autocovSurface,matindex)
#' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
#' grid()
#' }
.estimate_iid_distr_Imhof <- function(X, autocovs, l2norms){
  X <- .check_data(X)

  # # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
  # # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
  # means <- matrix(rep(rowMeans(X$data),ncol(X$data)), nrow=nrow(X$data))
  # l <- .obtain_autocov_eigenvalues(X$data - means,X$intraobs)
  l <- .obtain_autocov_eigenvalues(center(X))

  nl <- length(l)
  Reig <- .estimate_iid_distr_MC(X, autocovs, l2norms)$Reig
  ex <- seq(from = 0, to = max(Reig), length.out = 250)
  x <- ex * ncol(X$data)
  nx <- length(x)

  # Compute products of eigenvalues
  L <- rep(0,nl*nl)
  k <- 1
  for (i in l) {
    for (j in l){
      L[k] <- i*j
      k <- k+1
    }
  }

  # Obtain ecdf using Imhof function
  cdf <- rep(0,nx)
  for (i in 1:nx){
    cdf[i] <- 1 - CompQuadForm::imhof(x[i], L)$Qq
  }

  # ex <- w
  # ef <- t(cdf)
  # if(figure){
  #   # Check if any additional plotting parameters are present
  #   arguments <- list(...)
  #   if("main" %in% names(arguments)){
  #     graphics::plot(ex,ef,type = "l",...)
  #   }else{
  #     graphics::plot(ex,ef,type = "l",main = "ecdf obtained by Imhof estimation")
  #   }
  #   graphics::grid()
  # }

  list('ex'=ex, 'ef'=t(cdf))
}


#' Plot the autocorrelation function of a given FTS
#'
#' Plot a visual representation of the autocorrelation function of a given
#'   functional time series, including the upper i.i.d. bound.
#'
#' @param rho Autocorrelation values for each lag of
#' the functional time series obtained by calling the
#' function \code{obtain_FACF}.
#' @param SWN The upper prediction bound for the strong white noise iid
#'  distribution obtained by calling the function \code{.compute_FACF}.
#' @param WWN The upper prediction bound for the weak white noise iid
#'  distribution obtained by calling the function \code{.compute_FACF}.
#' @param ... Further arguments passed to the \code{plot} function.
#'
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' bbridge <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 15
#' upper_bound <- 0.95
#' fACF <- .compute_FACF(X = bbridge, lag.max = lag.max,
#'                       alpha=upper_bound, figure = FALSE)
#' .plot_FACF(rho = fACF$acfs,SWN = fACF$SWN_bound, WWN = fACF$WWN_bound)
#'
#' \donttest{
#' # Example 2
#'
#' N <- 200
#' v <- seq(from = 0, to = 1, length.out = 30)
#' sig <- 2
#' bbridge <- generate_brownian_bridge(N, v, sig)
#' lag.max <- 15
#' upper_bound <- 0.95
#' fACF <- .compute_FACF(X = bbridge, lag.max = lag.max,
#'                       alpha = upper_bound, figure = FALSE)
#' .plot_FACF(rho = fACF$acfs,SWN = fACF$SWN_bound,WWN = fACF$WWN_bound)
#' }
.plot_FACF <- function(rho, SWN, WWN, ...){
  # Define suitable lwd for plotting
  lag.max <- length(rho)

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))  arguments$xlab <- "Lag"
  if(!"ylab" %in% names(arguments))  arguments$ylab <- "ACF"
  if(!"ylim" %in% names(arguments))  arguments$ylim <- c(0, min(max(rho)*1.5,1))
  if(!"lwd"  %in% names(arguments))   arguments$lwd <- 1#lwd_1
  if(!"las"  %in% names(arguments))   arguments$las <- 1
  if(!"lend" %in% names(arguments))  arguments$lend <- 2
  if(!"yaxs" %in% names(arguments))  arguments$yaxs <- "i"
  if(!"xaxs" %in% names(arguments))  arguments$xaxs <- "i"
  if(!"main" %in% names(arguments))  arguments$main <- ""
  if(!"xlim" %in% names(arguments))  arguments$xlim <- c(0, length(rho)+1)
  arguments$x <- seq(1, length(rho), by = 1)
  arguments$y <- rho
  arguments$type <- "h"
  # arguments$xlim <- c(0,1)
  #arguments$ylim <- c(0,1.0)

  do.call(graphics::plot, arguments)
  #grid(lty = 1)
  do.call(graphics::lines, arguments)
  graphics::lines(x = arguments$x,
                  y = arguments$y,
                  type = arguments$type,
                  col = 'black',#"lightgrey",
                  # lwd = arguments$lwd - 2,
                  lend = 2)
  blue_col <- "#0073C2FF"
  graphics::abline(h = SWN, col = blue_col, lty = 2)# lwd = 4, lty = 2)
  graphics::lines(x = arguments$x, y = WWN, col = 'red', lty = 2)# lwd = 4, lty = 2)
  # graphics::legend(
  #   x = "topleft",
  #   legend = c(paste("i.i.d. bound (",ci*100," % conf.)",sep="")),
  #   col = blue_col,
  #   lty = 2,
  #   lwd = 4)
  box()
}


#' Obtain the partial autocorrelation function for a given FTS.
#'
#' Estimate the partial autocorrelation function for a given functional
#'  time series and its distribution under the hypothesis of strong
#'  functional white noise.
#'
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Discretization points of the curves.
#' @param lag.max Number of lagged covariance operators
#' of the functional time series that will be used
#' to estimate the partial autocorrelation function.
#' @param n_pcs Number of principal components
#' that will be used to fit the ARH(p) models.
#' @param ci A value between 0 and 1 that indicates
#' the confidence interval for the i.i.d. bounds
#' of the partial autocorrelation function. By default
#' \code{ci = 0.95}.
#' @param estimation Character specifying the
#' method to be used when estimating the distribution
#' under the hypothesis of functional white noise.
#' Accepted values are:
#' \itemize{
#'    \item "MC": Monte-Carlo estimation.
#'    \item "Imhof": Estimation using Imhof's method.
#' }
#' By default, \code{estimation = "MC"}.
#' @param figure Logical. If \code{TRUE}, plots the
#' estimated partial autocorrelation function with the
#' specified i.i.d. bound.
#' @param ... Further arguments passed to the \code{plot_FACF}
#' function.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{Blueline}: The upper prediction
#'     bound for the i.i.d. distribution.
#'     \item \code{rho}: Partial autocorrelation
#'     coefficients for
#'     each lag of the functional time series.
#' }
#' @export
#'
#' @references
#' Mestre G., Portela J., Rice G., Muñoz San Roque A., Alonso E. (2021).
#' \emph{Functional time series model identification and diagnosis by
#' means of auto- and partial autocorrelation analysis.}
#' Computational Statistics & Data Analysis, 155, 107108.
#' \url{https://doi.org/10.1016/j.csda.2020.107108}
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 5)
#' sig <- 2
#' set.seed(15)
#' X <- generate_brownian_bridge(N, v, sig)
#' .compute_FPACF(X,lag.max = 10, n_pcs = 2)
#'
#' \donttest{
#' # Example 2
#'
#' .compute_FPACF(X = funts(electricity), lag.max = 30, n_pcs = 5,
#'                alpha = 0.90, figure = TRUE)
#' }
.compute_FPACF <- function(X, n_pcs = NULL, lag.max = NULL,
                           alpha=0.95, figure = TRUE, ...){
  X <- .check_data(X)

  res <- nrow(X$data) #dv <- length(v)
  nobs <- ncol(X$data) #dt <- nrow(y)

  if(is.null(n_pcs)){
    max_pc <- 20

    # If there are less discretization points than max_pc, use the disc points
    num_fpc <- min(c(length(X$intraobs), max_pc))

    pca <- stats::princomp(t(X$data))#$scores[,1:num_fpc]
    eigs <- pca$sdev^2
    varprop <- as.numeric(cumsum(eigs[1:num_fpc]) / sum(eigs))

    # Select 95% Coverage or warn
    if(any(varprop>=0.95)){
      n_pcs <- which(varprop>=0.95)[1]
    }else{
      warning(paste0("Using ", num_fpc,
                     " functional principal components only explains ",
                     format(varprop[length(varprop)],digits = 3),"%",
                     " of the variance"))
      n_pcs <- num_fpc
    }
  }

  if( is.null(lag.max) ) lag.max <- 20

  # Initialize FPACF vector
  FPACF <- rep(NA, lag.max)

  FACF <- .compute_FACF(
    X = X, lag.max = 1,
    alpha = alpha, figure = FALSE, WWN = FALSE, ...)

  FPACF[1] <- FACF$acfs[1]

  vector_PACF <- FPACF

  show_varprop <- TRUE

  # Start loop for fitting ARH(p-1)
  for(pp in 2:lag.max){
    lag_PACF <- pp

    # 1 - Fit ARH(1) to the series
    if(show_varprop){
      Xest_ARIMA <-
        .fit_ARHp_FPCA(X=X, p = lag_PACF-1,
                       n_pcs = n_pcs)$x_est
      show_varprop <- FALSE
    }else{
      Xest_ARIMA <-
        .fit_ARHp_FPCA(X=X, p = lag_PACF-1,
                       n_pcs = n_pcs, show_varprop = FALSE)$x_est
    }

    # 2 - Fit ARH(1) to the REVERSED series
    x_rev <- X
    x_rev$data <- X$data[,seq(from = ncol(X$data), to = 1, by = -1)]

    Xest_ARIMA_REV <- .fit_ARHp_FPCA(
      X = x_rev, p = lag_PACF-1,
      n_pcs = n_pcs, show_varprop = FALSE)$x_est

    # 3 - Estimate covariance surface for PACF
    Xest_1 <- Xest_ARIMA
    Xest_2 <-
      Xest_ARIMA_REV[,seq(from = ncol(X$data), to = 1, by = -1)]

    res_filt_1 <- X$data - Xest_1
    res_filt_2 <- X$data - Xest_2

    # Cross-covariance surface
    sup_cov <- matrix(0, res, res)
    ini_serie <- max(which(is.na(Xest_1[1,]))) + 2
    fin_serie <- min(which(is.na(Xest_2[1,]))) - 2

    count <- 0
    for(jj in ini_serie:fin_serie){
      epsilon_1 <- as.matrix(res_filt_1[,jj])
      epsilon_2 <- as.matrix(res_filt_2[,jj - lag_PACF])

      sup_cov <- sup_cov + ( epsilon_1 %*% t(epsilon_2) )

      count <- count + 1
    }
    sup_cov <- sup_cov / count

    # Estimate traces
    var_1 <- matrix(0, res, res)
    count <- 0
    for (jj in ini_serie:nobs){
      epsilon_1 <- as.matrix(res_filt_1[, jj])

      var_1 <- var_1 + ( epsilon_1 %*% t(epsilon_1) )

      count <- count + 1
    }
    var_1 <- var_1 / count

    traza_1 <- dot_integrate_uneven(r = X$intraobs, v = diag(var_1))

    var_2 <- matrix(0, res, res)
    count <- 0
    for (jj in 1:fin_serie){
      epsilon_1 <- as.matrix(res_filt_2[,jj])
      epsilon_2 <- as.matrix(res_filt_2[,jj])

      var_2 <- var_2 + ( epsilon_1 %*% t(epsilon_2) )

      count <- count + 1
    }
    var_2 <- var_2 / count

    traza_2 <- dot_integrate_uneven(r = X$intraobs, v = diag(var_2))

    sup_corr <- sup_cov / ( sqrt(traza_1)*sqrt(traza_2) )

    vector_PACF[lag_PACF] <-
      sqrt( .obtain_suface_L2_norm(X$intraobs, list(Lag0 = sup_corr)) )
  }

  if(figure){
    .plot_FACF(rho = vector_PACF, SWN = FACF$SWN_bound, WWN=NULL, ylab='PACF', ...)
  }

  list(pacfs = vector_PACF, SWN = FACF$SWN_bound)
}


#' Fit an ARH(p) to a given functional time series
#'
#' Fit an \eqn{ARH(p)} model to a given functional time series. The fitted
#'   model is based on the model proposed in (Aue et al, 2015), first
#'   decomposing the original functional observations into a vector time series
#'   of \code{n_pcs} FPCA scores, and then fitting a vector autoregressive
#'   model of order \eqn{p} (\eqn{VAR(p)}) to the time series of the scores.
#'   Once fitted, the Karhunen-Loève expansion is used to re-transform the
#'   fitted values into functional observations.
#'
#' @inheritParams .compute_FPACF
#' @param p Numeric value specifying the order
#' of the functional autoregressive
#' model to be fitted.
#' @param show_varprop Logical. If \code{show_varprop = TRUE}, a plot of
#'   the proportion of variance explained by the first \code{n_pcs} functional
#'   principal components will be shown. By default \code{show_varprop = TRUE}.
#'
#' @export
#'
#' @examples
#' # Example 1
#'
# Simulate an ARH(1) process
#' N <- 250
#' dv <- 20
#' v <- seq(from = 0, to = 1, length.out = 20)
#'
#' phi <- 1.3 * ((v) %*% t(v))
#'
#' persp(v,v,phi,
#'       ticktype = "detailed",
#'       main = "Integral operator")
#'
#' set.seed(3)
#' white_noise <-  generate_brownian_bridge(N, v = v)
#'
#' y <- matrix(nrow = dv, ncol = N)
#' y[,1] <- white_noise$data[,1]
#' for(jj in 2:N){
#'     y[,jj] <- white_noise$data[,jj];
#'
#'     y[,jj] <- y[,jj] + .integral_operator(operator_kernel = phi,
#'                                     v = v, curve = y[,jj-1])
#' }
#'
#' # Fit an ARH(1) model
#' mod <- .fit_ARHp_FPCA(X = funts(y,intraobs = v), p = 1, n_pcs = 5)
#'
#' # Plot results
#' plot(v, y[,50], type = "l", lty = 1, ylab = "")
#' lines(v, mod$x_est[,50], col = "red")
#' legend("bottomleft", legend = c("real","est"),
#'        lty = 1, col = c(1,2))
#'
#' @references Aue, A., Norinho, D. D., Hormann, S. (2015).
#'  \emph{On the Prediction of Stationary Functional Time Series}
#'  Journal of the American Statistical Association,
#'  110, 378--392. \url{https://doi.org/10.1080/01621459.2014.909317}
.fit_ARHp_FPCA <- function(X, p, n_pcs, show_varprop = TRUE){
  X <- .check_data(X)

  nobs <- ncol(X$data) #dt <- nrow(y)
  res <- nrow(X$data) #dv <- length(v)

  # Step 1: FPCA decomposition of the curves
  pca <- stats::princomp(t(X$data))
  eigs <- pca$sdev^2
  varprop <- as.numeric(cumsum(eigs[1:n_pcs]) / sum(eigs))

  if(show_varprop){
    graphics::plot(1:n_pcs, varprop,
                   type = "b",
                   pch = 20,
                   main = "% Variance explained by FPCA",
                   xlab = "Number of components",
                   ylab = "% Var. Expl")
  }

  x_scores <- pca$scores[,1:n_pcs]


  # Step 2: Fit a VAR(p) model to the scores series
  x_scores_aux <- as.data.frame(x_scores)
  names(x_scores_aux) <- paste0("x",1:ncol(x_scores))
  mod <- vars::VAR(x_scores_aux, p = p)

  if(show_varprop){
    summary(mod)
  }

  fitted_vals <- stats::fitted(mod)

  # Fill with NA
  fitted_vals <- rbind(
    matrix(NA, nrow = nrow(x_scores)-nrow(fitted_vals),
           ncol = ncol(x_scores)),
    fitted_vals)

  # Step 3: reconstruct the curves
  x_rec <- matrix(nrow = res, ncol = nobs)
  for(ii in 1:nobs){
    x_rec[,ii] <- pca$center + pca$loadings[,1:n_pcs] %*% fitted_vals[ii,]
  }

  list(x_est = x_rec, mod = mod, fpca = pca,
       fitted_vals = fitted_vals, X = X)
}


#' Integral transformation of a curve using an integral operator
#'
#' Compute the integral transform of the curve \eqn{Y_i} with respect to
#'  a given integral operator \eqn{\Psi}. The transformation is given by
#'  \deqn{\Psi(Y_{i})(v) = \int \psi(u,v)Y_{i}(u)du}
#'
#' @param operator_kernel Matrix with the values of the kernel surface of
#'  the integral operator. The dimension of the matrix is \eqn{(g x m)},
#'  where \eqn{g} is the number of discretization points of the input curve
#'  and \eqn{m} is the number of discretization points of the output curve.
#' @param curve Vector containing the discretized values of a functional
#'  observation. The dimension of the matrix is \eqn{(1 x m)}, where
#'  \eqn{m} is the number of points observed in the curve.
#' @param v Numerical vector specifying the discretization points of the curves.
#' @return Returns a matrix the same size as \code{curve} with the transformed
#'  values.
#' @export
#'
#' @examples
#' # Example 1
#'
#' v <- seq(from = 0, to = 1, length.out = 20)
#' set.seed(10)
#' curve <- sin(v) + rnorm(length(v))
#' operator_kernel <- 0.6*(v %*% t(v))
#' hat_curve <- integral_operator(operator_kernel,curve,v)
integral_operator <- function(operator_kernel, curve, v){

  # Initialize output
  yhat <- rep(0, times = ncol(operator_kernel))

  # Perform the integral transformation
  for(ind_v in 1:ncol(operator_kernel)){
    yhat[ind_v] <- pracma::trapz(y = operator_kernel[,ind_v]*curve, x = v)
  }

  yhat
}
