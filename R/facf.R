#' Obtain the autocorrelation function for a given functional time series.
#'
#' @description Estimate the lagged autocorrelation function for a given
#' functional time series and its distribution under the
#' hypothesis of strong functional white noise. This graphic tool
#' can be used to identify seasonal patterns in the functional
#' data as well as auto-regressive or moving average terms.
#' i.i.d. bounds are included to test the presence of serial
#' correlation in the data.
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
#' @param ci A value between 0 and 1 that indicates
#' the confidence interval for the i.i.d. bounds
#' of the autocorrelation function. By default
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
#' estimated autocorrelation function with the
#' specified i.i.d. bound.
#' @param ... Further arguments passed to the \code{plot_FACF}
#' function.
#' @return Return a list with:
#' \itemize{
#'     \item \code{Blueline}: The upper prediction
#'     bound for the i.i.d. distribution.
#'     \item \code{rho}: Autocorrelation values for
#'     each lag of the functional time series.
#' }
#' @export
#' @references
#' Mestre G., Portela J., Rice G., Muñoz San Roque A., Alonso E. (2021).
#' \emph{Functional time series model identification and diagnosis by
#' means of auto- and partial autocorrelation analysis.}
#' Computational Statistics & Data Analysis, 155, 107108.
#' \url{https://doi.org/10.1016/j.csda.2020.107108}
#'
#' Mestre, G., Portela, J., Muñoz-San Roque, A., Alonso, E. (2020).
#' \emph{Forecasting hourly supply curves in the Italian Day-Ahead
#' electricity market with a double-seasonal SARMAHX model.}
#' International Journal of Electrical Power & Energy Systems,
#' 121, 106083. \url{https://doi.org/10.1016/j.ijepes.2020.106083}
#'
#' Kokoszka, P., Rice, G., Shang, H.L. (2017).
#' \emph{Inference for the autocovariance of a functional
#'  time series under conditional heteroscedasticity}
#' Journal of Multivariate Analysis,
#' 162, 32--50. \url{https://doi.org/10.1016/j.jmva.2017.08.004}
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 5)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' obtain_FACF(Y,v,20)
#'
#' \donttest{
#' # Example 2
#'
#' data(elec_prices)
#' v <- seq(from = 1, to = 24)
#' nlags <- 30
#' obtain_FACF(Y = as.matrix(elec_prices),
#' v = v,
#' nlags = nlags,
#' ci = 0.95,
#' figure = TRUE)
#' }
.compute_FACF <- function(X, v, lag.max = NULL, ci=0.95,
                          mc_estimation = TRUE,figure = TRUE,...){

  # v <- x$
  # X <- x$data

  if(is.null(lag.max))
    lag.max <- 20 # 10 * log(ncol(X)/, base=10)

  # Autocov surfaces
  #autocovSurface <- obtain_autocovariance(Y,nlags)
  autocovs <- .compute_autocovariance(X, lag.max)

  # L2 norm autocov surfaces
  #matindex <- obtain_suface_L2_norm (v,autocovSurface)
  #matindex <- matindex[-1]
  l2norms <- .obtain_suface_L2_norm(v, autocovs)
  l2norms <- l2norms[-1] # Drop Lag 0

  # Estimate distribution of iid bound
  if(mc_estimation){
    # set.seed(123)
    # iid.distribution <-
    #   estimate_iid_distr_MC(Y,v,autocovSurface,matindex,
    #                         figure = FALSE)
    # set.seed(123)
    iid.distribution <-
      .estimate_iid_distr_MC(X, v, autocovs, l2norms)
  } else{
    # set.seed(123)
    # iid.distribution <-
    #   estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex,figure = FALSE)
    #
    # set.seed(123)
    iid.distribution <-
      .estimate_iid_distr_Imhof(X, v, autocovs, l2norms)
  }

  # Obtain autocorrelation
  normalization.value <-
    pracma::trapz(v, diag(autocovs$Lag0))
  rho <- sqrt(l2norms) / normalization.value

  # Obtain iid bound for specified confidence value
  cutoff <- rep(NA,length(ci))
  for(ii in 1:length(ci)){
    idx <- min(which(iid.distribution$ef >= ci[ii]))
    cutoff[ii] <-
      sqrt( iid.distribution$ex[ idx ] ) /
      normalization.value
  }

  if(figure){
    .plot_FACF(rho,cutoff, ci, ...)
  }

  list('cutoff'=cutoff, 'rho'=rho)
}


#' Estimate Autocovariance Function of funts
#'
#' Obtain the empirical autocovariance function for
#' lags \eqn{= 0,...,}\code{nlags} of the functional time
#' series. Given \eqn{X_{1},...,X_{T}} a functional time
#' series, the sample autocovariance functions
#' \eqn{\hat{C}_{h}(u,v)} are given by:
#' \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
#' where
#' \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
#' denotes the sample mean function.
#'
#' @param X funts object
#' @param nlags Number of lagged covariance operators
#' of the functional time series that will be used
#' to estimate the autocorrelation function.
#'
#' @return Return a list with the lagged autocovariance
#' functions estimated from the data. Each function is given
#' by a \eqn{(m x m)} matrix, where \eqn{m} is the
#' resolution observed in each curve.
#' @export
#'
#' @examples
#' .compute_autocovariance(funts(electricity), 10)
.compute_autocovariance <- function(X, nlags){
  X_demean <- X - rowMeans(X) # TODO:: Check row means
  res <- nrow(X)
  obs <- ncol(X)

  autocovs <- list()
  for(nlag in 0:nlags){
    autocovs[[paste0("Lag", nlag)]] <-
      matrix(0, nrow = res, ncol = res)
    for( ind_curve in (1+nlag):obs ){
      autocovs[[paste0("Lag", nlag)]] <-
        autocovs[[paste0("Lag", nlag)]] +
        X_demean[,ind_curve - nlag] %*% t(X_demean[,ind_curve])
    }
    autocovs[[paste0("Lag", nlag)]] <-
      autocovs[[paste0("Lag", nlag)]] / (obs-1)
  }

  autocovs
}


#' Obtain L2 norm of the autocovariance functions
#'
#' Returns the L2 norm of the lagged autocovariance
#' functions \eqn{\hat{C}_{h}}. The L2 norm of these
#' functions is defined as
#' \deqn{\sqrt(\int \int \hat{C}^{2}_{h}(u,v)du dv)}
#'
#' @param v Discretization points of the curves, by default
#' \code{seq(from = 0, to = 1, length.out = 100)}.
#' @param autocovSurface An \eqn{(m x m)} matrix with
#' the discretized values of the autocovariance operator
#' \eqn{\hat{C}_{0}}, obtained by calling the
#' function \code{obtain_autocovariance}. The value
#' \eqn{m} indicates the number of points observed in
#' each curve.
#'
#' @return A vector containing the L2 norm of the
#' lagged autocovariance functions \code{autocovSurface}.
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 1
#' autocovSurface <- obtain_autocovariance(Y=Y,nlags = nlags)
#' norms <- obtain_suface_L2_norm(v = v,autocovSurface = autocovSurface)
#' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 1)
#' title(sub = paste0("Lag ",1," - L2 Norm: ",norms[2]))
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 2
#' autocovSurface <- obtain_autocovariance(Y=Y,nlags = nlags)
#' norms <- obtain_suface_L2_norm(v = v,autocovSurface = autocovSurface)
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
  nlags <- length(autocovs)
  matindex <- rep(NA, nlags)
  res <- nrow(autocovs$Lag0)

  # L2 norm of squared autocov surfaces
  for(rr in 1:nlags){
    autocov2 <- autocovs[[paste("Lag",rr-1,sep="")]]^2

    norm.vec <- matrix(NA, nrow = res, ncol = 1)
    for(ii in 1:res){
      norm.vec[ii] <- pracma::trapz(intraobs, autocov2[,ii])
    }

    norm.aux <- pracma::trapz(intraobs, norm.vec)
    matindex[rr] <- norm.aux
  }

  matindex
}


#' Estimate distribution of the fACF under the iid. hypothesis using MC method
#'
#' Estimate the distribution of the autocorrelation function
#' under the hypothesis of strong functional white noise. This
#' function uses Montecarlo's method to estimate the distribution.
#'
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Discretization points of the curves, by default
#' \code{seq(from = 0, to = 1, length.out = 100)}.
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
#' @param figure Logical. If \code{TRUE}, plots the
#' estimated distribution.
#' @param ... Further arguments passed to the  \code{plot}
#' function.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{ex}: Knots where the
#'     distribution has been estimated
#'     \item \code{ef}: Discretized values of
#'     the estimated distribution.
#'     \item \code{Reig}: Raw values of the i.i.d.
#'     statistic for each MC simulation.
#' }
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 1
#' autocovSurface <- obtain_autocovariance(Y,nlags)
#' matindex <- obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' MC_dist <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex)
#' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
#' grid()
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 20
#' autocovSurface <- obtain_autocovariance(Y,nlags)
#' matindex <- obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' MC_dist <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex)
#' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
#' grid()
#' }
.estimate_iid_distr_MC <-
  function(X, v, autocovSurface, matindex, nsims= 10000){

    # TODO:: Update Means
    # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
    # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
    means <- matrix(rep(rowMeans(X),ncol(X)), nrow=nrow(X))
    l <- .obtain_autocov_eigenvalues(X - means,v)

    neig <- length(l)
    Reig <- rep(0,nsims)

    for(jj in 1:neig){
      for(kk in 1:neig){
        Reig <- Reig +
          stats::rchisq(n = nsims,df = 1)*l[jj]*l[kk]
      }
    }
    #Reig=Reig/nrow(Y)
    Reig <- Reig/ncol(X)

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
#' function \eqn{\hat{C}_{0}}. This functions returns the
#' eigenvalues which are greater than the value \code{epsilon}.
#'
#' @param v Discretization points of the curves, by default
#' \code{seq(from = 0, to = 1, length.out = 100)}.
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param epsilon Value used to determine how many
#' eigenvalues will be returned. The eigenvalues
#' \eqn{\lambda_{j}/\lambda_{1} > \code{epsilon}}
#' will be returned.
#' By default \code{epsilon = 0.0001}.
#'
#' @return A vector containing the \eqn{k} eigenvalues
#' greater than \code{epsilon}.
#' @export
#'
#' @examples
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' lambda <- obtain_autocov_eigenvalues(v = v, Y = Y)
.obtain_autocov_eigenvalues <- function(X, v, epsilon = 0.0001){

  nobs <- ncol(X) # nt
  res <- nrow(X) # nv

  # w(i,j) = integral of the product of the curves i and j
  W <- matrix(0,nrow = nobs,ncol = nobs)

  for(ii in 1:nobs){
    mat.aux <- matrix(rep(X[,ii],each = nobs),
                      nrow = nobs, ncol = res)*t(X)
    for(jj in 1:nobs){
      W[ii,jj] = pracma::trapz(v,mat.aux[jj,])
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
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Discretization points of the curves, by default
#' \code{seq(from = 0, to = 1, length.out = 100)}.
#' @param autocovSurface An \eqn{(m x m)} matrix with the discretized
#' values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained
#' by calling the function \code{obtain_autocovariance}.
#' The value \eqn{m} indicates the number of points observed
#' in each curve.
#' @param matindex A vector containing the L2 norm of
#' the autocovariance function. It can be obtained by calling
#' function \code{obtain_suface_L2_norm}.
#' @param figure Logical. If \code{TRUE}, plots the
#' estimated distribution.
#' @param ... Further arguments passed to the  \code{plot}
#' function.
#'
#' @return Return a list with:
#' \itemize{
#'     \item \code{ex}: Knots where the
#'     distribution has been estimated
#'     \item \code{ef}: Discretized values of
#'     the estimated distribution.
#' }
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 1
#' autocovSurface <- obtain_autocovariance(Y,nlags)
#' matindex <- obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' Imhof_dist <- estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex)
#' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
#' grid()
#'
#' \donttest{
#' # Example 2
#'
#' N <- 400
#' v <- seq(from = 0, to = 1, length.out = 50)
#' sig <- 2
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' autocovSurface <- obtain_autocovariance(Y,nlags)
#' matindex <- obtain_suface_L2_norm (v,autocovSurface)
#' # Remove lag 0
#' matindex <- matindex[-1]
#' Imhof_dist <- estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex)
#' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
#' grid()
#' }
.estimate_iid_distr_Imhof <- function(X,v,autocovs,l2norms){


  # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
  # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
  means <- matrix(rep(rowMeans(X),ncol(X)), nrow=nrow(X))
  l <- .obtain_autocov_eigenvalues(X - means,v)

  nl <- length(l)
  Reig <- .estimate_iid_distr_MC(X,v,autocovs,l2norms)$Reig
  ex <- seq(from = 0, to = max(Reig), length.out = 250)
  x <- ex * ncol(X)
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
#' @description Plot a visual representation of the autocorrelation function
#' of a given functional time series, including the upper i.i.d.
#' bound.
#'
#' @param rho Autocorrelation values for each lag of
#' the functional time series obtained by calling the
#' function \code{obtain_FACF}.
#' @param Blueline The upper prediction bound for the
#' i.i.d. distribution obtained by calling the
#' function \code{obtain_FACF}.
#' @param ci Value between 0 and 1 that was used
#' when calling the function \code{obtain_FACF}.
#' This value is only used to display information
#' in the figure.
#' @param ... Further arguments passed to the  \code{plot}
#' function.
#'
#' @export
#'
#' @examples
#' # Example 1
#'
#' N <- 100
#' v <- seq(from = 0, to = 1, length.out = 10)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 15
#' upper_bound <- 0.95
#' fACF <- obtain_FACF(Y = bbridge,v = v,nlags = nlags,ci=upper_bound,figure = FALSE)
#' plot_FACF(rho = fACF$rho,Blueline = fACF$Blueline,ci = upper_bound)
#'
#' \donttest{
#' # Example 2
#'
#' N <- 200
#' v <- seq(from = 0, to = 1, length.out = 30)
#' sig <- 2
#' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
#' nlags <- 15
#' upper_bound <- 0.95
#' fACF <- obtain_FACF(Y = bbridge,v = v,nlags = nlags,ci=upper_bound,figure = FALSE)
#' plot_FACF(rho = fACF$rho,Blueline = fACF$Blueline,ci = upper_bound)
#' }
.plot_FACF <- function(rho, cutoff, ci, ...){
  # Define suitable lwd for plotting
  nlags <- length(rho)
  # if(nlags <= 30){
  #   lwd_1 <- 8
  # }else if(nlags <= 50){
  #   lwd_1 <- 6
  # }else if(nlags <= 150){
  #   lwd_1 <- 4
  # }else{
  #   lwd_1 <- 3
  # }

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))  arguments$xlab <- "Lag"
  if(!"ylab" %in% names(arguments))  arguments$ylab <- "ACF"
  if(!"ylim" %in% names(arguments))  arguments$ylim <- c(0, max(rho)*1.5)
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
  graphics::abline(h = cutoff, col = blue_col, lty = 2)# lwd = 4, lty = 2)
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
#' @description Estimate the partial autocorrelation
#' function for a given functional time series and its
#' distribution under the hypothesis of strong functional
#' white noise.
#'
#' @param Y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Discretization points of the curves.
#' @param nlags Number of lagged covariance operators
#' of the functional time series that will be used
#' to estimate the partial autocorrelation function.
#' @param n_harm Number of principal components
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
#' Y <- simulate_iid_brownian_bridge(N, v, sig)
#' obtain_FPACF(Y,v,10, n_harm = 2)
#'
#' \donttest{
#' # Example 2
#'
#' data(elec_prices)
#' v <- seq(from = 1, to = 24)
#' nlags <- 30
#' obtain_FPACF(Y = as.matrix(elec_prices),
#' v = v,
#' nlags = nlags,
#' n_harm = 5,
#' ci = 0.95,
#' figure = TRUE)
#' }
.compute_FPACF <- function(X, v, n_harm = NULL, lag.max = NULL,
                           ci=0.95, figure = TRUE, ...){

  x <- X
  res <- nrow(x) #dv <- length(v)
  nobs <- ncol(x) #dt <- nrow(y)

  if(is.null(n_harm)){
    max_pc <- 20

    # If there are less discretization points than max_pc, use the disc points
    num_fpc <- min(c(length(v), max_pc))

    x_fd <- mat2fd(mat_obj = t(X),argvals = v, range_val = range(v))
    pca <- fda::pca.fd(fdobj = x_fd, nharm = num_fpc)

    varprop <- cumsum(pca$varprop)

    # Select 95% Coverage or warn
    if(any(varprop>=0.95)){
      n_harm <- which(varprop>=0.95)[1]
    }else{
      warning(paste0("Using ",
                     num_fpc,
                     " functional principal components only explains ",
                     format(varprop[length(varprop)],digits = 3),"%",
                     " of the variance"))
      n_harm <- num_fpc
    }
  }

  if(is.null(lag.max))
    lag.max <- 20

  # Initialize FPACF vector
  FPACF <- rep(NA, lag.max)

  FACF <- .compute_FACF(
    X = x, v = v, lag.max = 1,
    ci = ci, figure = FALSE, ...)

  cutoff <- FACF$cutoff
  FPACF[1] <- FACF$rho[1]

  vector_PACF <- FPACF

  show_varprop = T

  # Start loop for fitting ARH(p-1)
  for(pp in 2:lag.max){
    lag_PACF <- pp

    # 1 - Fit ARH(1) to the series
    if(show_varprop){
      Xest_ARIMA <-
        .fit_ARHp_FPCA(x = x, v = v, p = lag_PACF-1,
                       n_harm = n_harm)$x_est
      show_varprop <- FALSE
    }else{
      Xest_ARIMA <-
        .fit_ARHp_FPCA(x = x, v = v, p = lag_PACF-1,
                       n_harm = n_harm, show_varprop = F)$x_est
    }

    # if(F){
    #   kkk <- 30
    #   graphics::plot(v,y[kkk,],type = "b")
    #   graphics::lines(v,Yest_ARIMA[kkk,],col = "red")
    # }

    # 2 - Fit ARH(1) to the REVERSED series
    x_rev <- x[,seq(from = ncol(x), to = 1, by = -1)]

    Xest_ARIMA_REV <- .fit_ARHp_FPCA(
      x = x_rev, v = v, p = lag_PACF-1,
      n_harm = n_harm, show_varprop = F)$x_est

    # 3 - Estimate covariance surface for PACF
    Xest_1 <- Xest_ARIMA
    Xest_2 <-
      Xest_ARIMA_REV[,seq(from = ncol(x), to = 1, by = -1)]

    res_filt_1 = x - Xest_1
    res_filt_2 = x - Xest_2

    # Cross-covariance surface
    sup_cov <- matrix(0, res, res)
    ini_serie <- max(which(is.na(Xest_1[1,])))+2
    fin_serie = min(which(is.na(Xest_2[1,])))-2

    count <- 0
    for(jj in ini_serie:fin_serie){
      epsilon_1 <- as.matrix(res_filt_1[,jj])
      epsilon_2 <- as.matrix(res_filt_2[,jj-lag_PACF])

      sup_cov <- sup_cov + ( epsilon_1 %*% t(epsilon_2) )

      count <- count + 1
    }

    sup_cov <- sup_cov / count

    # if(F){
    #   graphics::persp(v,v,sup_cov,theta = 330, phi = 20,
    #                   main = paste0("Covariance surface PACF lag ",lag_PACF),
    #                   ticktype = "detailed")
    # }

    # # L2 norm of covariance surface
    # if(F){
    #   sqrt(fdaACF::obtain_suface_L2_norm(v, list(Lag0 = sup_cov)))
    # }


    # Estimate traces
    var_1 <- matrix(0, res, res)
    count <- 0
    for (jj in ini_serie:nobs){
      epsilon_1 <- as.matrix(res_filt_1[,jj])

      var_1 <- var_1 + ( epsilon_1 %*% t(epsilon_1) )

      count <- count + 1
    }
    var_1 <- var_1 / count

    # if(F){
    #   graphics::persp(v,v,var_1,theta = 330, phi = 20, ticktype = "detailed", main = "var 1")
    # }
    traza_1 <- pracma::trapz(v, diag(var_1))

    var_2 <- matrix(0, res, res)
    count <- 0
    for (jj in 1:fin_serie){
      epsilon_1 <- as.matrix(res_filt_2[,jj])
      epsilon_2 <- as.matrix(res_filt_2[,jj])

      var_2 <- var_2 + ( epsilon_1 %*% t(epsilon_2) )

      count <- count + 1


    }
    var_2 <- var_2 / count
    # if(F){
    #   graphics::persp(v,v,var_2,theta = 330, phi = 20, ticktype = "detailed", main = "var 2")
    # }

    traza_2 <- pracma::trapz(v, diag(var_2))

    sup_corr <- sup_cov / ( sqrt(traza_1)*sqrt(traza_2) )

    # # L2 norm of cross-correlation surface
    # if(F){
    #   sqrt(fdaACF::obtain_suface_L2_norm(v, list(Lag0 = sup_corr)))
    # }

    vector_PACF[lag_PACF] <-
      sqrt( .obtain_suface_L2_norm(v, list(Lag0 = sup_corr)) )
  }

  if(figure){
    .plot_FACF(vector_PACF, cutoff, ci, ...)
  }

  list(cutoff = cutoff, rho = vector_PACF)
}


#' Fit an ARH(p) to a given functional time series
#'
#' @description Fit an \eqn{ARH(p)} model to a given functional
#' time series. The fitted model is based on the model proposed
#' in (Aue et al, 2015), first decomposing the original
#' functional observations into a vector time series of \code{n_harm}
#' FPCA scores, and then fitting a vector autoregressive
#' model of order \eqn{p} (\eqn{VAR(p)}) to the time series
#' of the scores. Once fitted, the Karhunen-Loève expansion
#' is used to re-transform the fitted values into functional
#' observations.
#'
#' @param y Matrix containing the discretized values
#' of the functional time series. The dimension of the
#' matrix is \eqn{(n x m)}, where \eqn{n} is the
#' number of curves and \eqn{m} is the number of points
#' observed in each curve.
#' @param v Numeric vector that contains the
#' discretization points of the curves.
#' @param p Numeric value specifying the order
#' of the functional autoregressive
#' model to be fitted.
#' @param n_harm Numeric value specifying the number
#' of functional principal components to be used when fitting
#' the \eqn{ARH(p)} model.
#' @param show_varprop Logical. If \code{show_varprop = TRUE},
#' a plot of the proportion of variance explained by the first
#' \code{n_harm} functional principal components will be shown.
#' By default \code{show_varprop = TRUE}.
#'
#' @export
#'
#' @examples
#' # Example 1
#'
#' # Simulate an ARH(1) process
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
#' white_noise <-  simulate_iid_brownian_bridge(N, v = v)
#'
#' y <- matrix(nrow = N, ncol = dv)
#' y[1,] <- white_noise[1,]
#' for(jj in 2:N){
#'     y[jj,] <- white_noise[jj,];
#'
#'     y[jj,]=y[jj,]+integral_operator(operator_kernel = phi,
#'                                     v = v,
#'                                     curve = y[jj-1,])
#' }
#'
#' # Fit an ARH(1) model
#' mod <- fit_ARHp_FPCA(y = y,
#'                      v = v,
#'                      p = 1,
#'                      n_harm = 5)
#'
#' # Plot results
#' plot(v, y[50,], type = "l", lty = 1, ylab = "")
#' lines(v, mod$y_est[50,], col = "red")
#' legend("bottomleft", legend = c("real","est"),
#'        lty = 1, col = c(1,2))
#'
#' @references Aue, A., Norinho, D. D., Hormann, S. (2015).
#' \emph{On the Prediction of Stationary Functional
#' Time Series}
#' Journal of the American Statistical Association,
#' 110, 378--392. \url{https://doi.org/10.1080/01621459.2014.909317}
.fit_ARHp_FPCA <- function(x, v, p, n_harm, show_varprop = T){
  nobs <- ncol(x) #dt <- nrow(y)
  res <- nrow(x) #dv <- length(v)

  # Step 1: FPCA decomposition of the curves
  x_fd <- mat2fd(mat_obj = t(x), argvals = v,
                 range_val = range(v))
  pca <- fda::pca.fd(fdobj = x_fd, nharm = n_harm)

  if(show_varprop){
    graphics::plot(1:n_harm, cumsum(pca$varprop),
                   type = "b",
                   pch = 20,
                   main = "% Variance explained by FPCA",
                   xlab = "Number of components",
                   ylab = "% Var. Expl")
  }

  x_scores <- pca$scores


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

  # if(F){
  #   graphics::plot(y_scores[,1], type = "l")
  #   graphics::lines(fitted_vals[,1], col = "red")
  # }

  # Step 3: reconstruct the curves
  x_rec <- matrix(nrow = res, ncol = nobs)
  for(ii in 1:nobs){
    fd_aux <- .reconstruct_fd_from_PCA(
      pca_struct = pca, scores = fitted_vals[ii,],
      centerfns = TRUE)

    fd_aux_mat <- fda::eval.fd(evalarg = v, fdobj = fd_aux)

    # if(F){
    #   graphics::plot(v,y[ii,], type = "b")
    #   graphics::lines(v,fd_aux_mat,col = "red")
    # }

    x_rec[,ii] <- fd_aux_mat
  }


  list(x_est = x_rec, mod = mod, fpca = pca,
       fitted_vals = fitted_vals, x = x, v = v)
}


#' Obtain the reconstructed curves after PCA
#'
#' @description This function reconstructs the
#' functional curves from a score vector
#' using the basis obtained after applying
#' functional PCA. This allows the user to
#' draw estimations from the joint density of
#' the FPCA scores and reconstruct the curves
#' for those new scores.
#'
#' @param pca_struct List obtained after calling
#' function \code{pca.fd}.
#' @param scores Numerical vector that contains the
#' scores of the fPCA decomposition for one
#' functional observation.
#' @param centerfns Logical value specifying
#' wheter the FPCA performed used \code{centerfns = T}
#' or \code{centerfns = F}. By default
#' \code{centerfns = T}.
#'
#' @return Returns a object of type \code{fd}
#' that contains the reconstructed curve.
#' @export
#'
#' @examples
#' # Example 1
#'
#' # Simulate fd
#' nobs <- 200
#' dv <- 10
#' basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=10)
#' set.seed(5)
#' C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
#' fd_sim <- fda::fd(coef=t(C),basis)
#'
#' # Perform FPCA
#' pca_struct <- fda::pca.fd(fd_sim,nharm = 6)
#'
#' # Reconstruct first curve
#' fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
#' plot(fd_sim[1])
#' plot(fd_rec, add = TRUE, col = "red")
#' legend("topright",
#'        legend = c("Real Curve", "PCA Reconstructed"),
#'        col = c("black","red"),
#'        lty = 1)
#'
#' # Example 2 (Perfect reconstruction)
#'
#' # Simulate fd
#' nobs <- 200
#' dv <- 7
#' basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=dv)
#' set.seed(5)
#' C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
#' fd_sim <- fda::fd(coef=t(C),basis)
#'
#' # Perform FPCA
#' pca_struct <- fda::pca.fd(fd_sim,nharm = dv)
#'
#' # Reconstruct first curve
#' fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
#' plot(fd_sim[1])
#' plot(fd_rec, add = TRUE, col = "red")
#' legend("topright",
#'        legend = c("Real Curve", "PCA Reconstructed"),
#'        col = c("black","red"),
#'        lty = 1)
.reconstruct_fd_from_PCA <- function(
    pca_struct, scores, centerfns = T){

  if (centerfns){
    out_curve <- pca_struct$meanfd
    for(rr in 1:ncol(pca_struct$scores)){
      out_curve <- out_curve + pca_struct$harmonics[rr]*scores[rr]
    }
  }else{
    out_curve <- pca_struct$harmonics[1]*scores[1]
    for(rr in 2:ncol(pca_struct$scores)){
      out_curve <- out_curve + pca_struct$harmonics[rr]*scores[rr]
    }
  }

  out_curve
}
