#' ACF/PACF Functions
#'
#' This function computes the ACF/PACF of data. This can be applied on traditional
#'  scalar time series or functional time series defined in [dfts()].
#'
#' @param x Object for computation of (partial) autocorrelation function
#'  (see \code{acf()} or \code{pacf}).
#' @param lag.max Number of lagged covariance estimators for the time series
#'  that will be used to estimate the (partial) autocorrelation function.
#' @param ... Additional parameters to appropriate function
#'
#' @seealso [stats::acf()]
#'
#' @name acf
#'
#' @return List with ACF or PACF values and plots
#' @export
#'
#' @examples
#' acf(1:10)
NULL


#' @rdname acf
#'
#' @export
acf <- function(x, lag.max = NULL, ...) UseMethod("acf")
#' @rdname acf
#'
#' @export
acf.default <- function(x, lag.max = NULL, ...) stats::acf(x)


#' @rdname acf
#'
#' @export
pacf <- function(x, lag.max = NULL, ...) UseMethod("pacf")
#' @rdname acf
#'
#' @export
pacf.default <- function(x, lag.max = NULL, ...) stats::pacf(x)


#' Obtain the autocorrelation function for a given functional time series.
#'
#' @param alpha A value between 0 and 1 that indicates significant level for
#'  the confidence interval for the i.i.d. bounds of the (partial) autocorrelation
#'  function. By default \code{alpha = 0.05}.
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
#' @param figure Logical. If \code{TRUE}, plots the estimated
#'  function with the specified bounds.
#' @param ... Further arguments passed to the \code{.plot_FACF} function.
#'
#' @return
#' \itemize{
#'     \item \code{acfs/pacfs}: Autocorrelation values for
#'     each lag of the functional time series.
#'     \item \code{SWN_bound}: The upper prediction
#'     bound for the i.i.d. distribution under strong white noise assumption.
#'     \item \code{WWN_bound}: The upper prediction
#'     bound for the i.i.d. distribution under weak white noise assumption.
#'     \item \code{plot}: Plot of autocorrelation values for
#'     each lag of the functional time series.
#' }
#'
#' @references Mestre G., Portela J., Rice G., Munoz San Roque A., Alonso E. (2021).
#'  \emph{Functional time series model identification and diagnosis by
#'  means of auto- and partial autocorrelation analysis.}
#'  Computational Statistics & Data Analysis, 155, 107108.
#'
#' @references Mestre, G., Portela, J., Munoz San Roque, A., Alonso, E. (2020).
#'  \emph{Forecasting hourly supply curves in the Italian Day-Ahead
#'  electricity market with a double-seasonal SARMAHX model.}
#'  International Journal of Electrical Power & Energy Systems,
#'  121, 106083.
#'
#' @references Kokoszka, P., Rice, G., Shang, H.L. (2017).
#'  \emph{Inference for the autocovariance of a functional
#'  time series under conditional heteroscedasticity}
#'  Journal of Multivariate Analysis,
#'  162, 32--50.
#'
#' @examples
#' x <- generate_brownian_bridge(100, seq(0,1,length.out=20))
#' acf(x,20)
#'
#' @export
#' @rdname acf
acf.dfts <- function(x, lag.max = NULL, alpha=0.05,
                      method = c('Welch','MC','Imhof'),
                      WWN = TRUE, figure = TRUE, ...){
  x <- dfts(x)

  if(is.null(lag.max))
    lag.max <- 20 # 10 * log(ncol(x)/, base=10)
  if(lag.max<1) stop('Increase lag.max to be greater than 0.',call. = FALSE)
  lag.max <- min(lag.max, ncol(x$data)-1)

  # Autocov surfaces
  # autocovs <- .compute_autocovariance(x, lag.max)
  autocovs <- autocovariance(x, 0:lag.max)

  # L2 norm autocov surfaces
  # l2norms <- .obtain_suface_L2_norm(x$fparam, autocovs)
  # l2norms <- l2norms[-1] # Drop Lag 0
  l2norms <- sapply(1:lag.max, function(idx,autocovs,res){
    dot_integrate(
      dot_integrate_col(v=t(autocovs[[idx+1]]^2),r=res),
      r=res)
  },autocovs=autocovs,res=x$fparam)

  # Obtain autocorrelation estimates
  normalization.value <-
    dot_integrate(r = x$fparam, v = diag(autocovs$Lag0))
  rho <- sqrt(l2norms) / normalization.value

  # Estimate distribution of SWN (iid) bound
  method <- .verify_input(method, c('welch','mc','imhof'))
  if(method=='welch'){
    # Obtain SWN bound for specified confidence value
    SWN_bound <- rep(NA,length(alpha))
    for(ii in 1:length(alpha)){
      SWN_bound[ii] <-
        sqrt(Q_WS_quantile_iid(x$data, alpha=alpha)$quantile) /
        ( sqrt(NCOL(x$data)) * normalization.value )
    }
  } else{
    if(method=='mc'){
      iid.distribution <- .estimate_iid_distr_MC(x, autocovs, l2norms)
    } else if(method=='imhof'){
      iid.distribution <- .estimate_iid_distr_Imhof(x, autocovs, l2norms)
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
    quantile <- Q_WS_quantile(x$data, lags, alpha=alpha, M=NULL, low_disc=FALSE)$quantile
    WWN_bound <- sqrt(quantile) /
      ( sqrt(ncol(x$data)) * normalization.value )
  }

  # Plot
  plt <- .plot_FACF(rho,SWN=SWN_bound, WWN=WWN_bound, ...)
  if(figure){
    plt
  }

  invisible( list('acfs'=rho, 'SWN_bound'=SWN_bound, 'WWN_bound'=WWN_bound, 'plot'=plt) )
}


#' Obtain the partial autocorrelation function for a given FTS.
#'
#' @param n_pcs Number of principal components that will be used to fit the
#'   ARH(p) models.
#'
#' @examples
#' x <- generate_brownian_bridge(100, seq(0,1,length.out=20))
#' pacf(x,lag.max = 10, n_pcs = 2)
#'
#' @export
#' @rdname acf
pacf.dfts <- function(x, lag.max = NULL, n_pcs = NULL,
                      alpha=0.95, figure = TRUE, ...){
  x <- dfts(x)

  res <- nrow(x$data) #dv <- length(v)
  nobs <- ncol(x$data) #dt <- nrow(y)

  # Allow TVE
  if(is.null(n_pcs)){
    max_pc <- 20

    # If there are less discretization points than max_pc, use the disc points
    num_fpc <- min(c(length(x$fparam), max_pc))

    pca <- stats::princomp(t(x$data))#$scores[,1:num_fpc]
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

  FACF <- acf.dfts(
    x = x, lag.max = 1,
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
        .fit_ARHp_FPCA(x=x, p = lag_PACF-1,
                       n_pcs = n_pcs)$x_est
      show_varprop <- FALSE
    }else{
      Xest_ARIMA <-
        .fit_ARHp_FPCA(x=x, p = lag_PACF-1,
                       n_pcs = n_pcs, show_varprop = FALSE)$x_est
    }

    # 2 - Fit ARH(1) to the REVERSED series
    x_rev <- x
    x_rev$data <- x$data[,seq(from = ncol(x$data), to = 1, by = -1)]

    Xest_ARIMA_REV <- .fit_ARHp_FPCA(
      x = x_rev, p = lag_PACF-1,
      n_pcs = n_pcs, show_varprop = FALSE)$x_est

    # 3 - Estimate covariance surface for PACF
    Xest_1 <- Xest_ARIMA
    Xest_2 <-
      Xest_ARIMA_REV[,seq(from = ncol(x$data), to = 1, by = -1)]

    res_filt_1 <- x$data - Xest_1
    res_filt_2 <- x$data - Xest_2

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

    traza_1 <- dot_integrate(r = x$fparam, v = diag(var_1))

    var_2 <- matrix(0, res, res)
    count <- 0
    for (jj in 1:fin_serie){
      epsilon_1 <- as.matrix(res_filt_2[,jj])
      epsilon_2 <- as.matrix(res_filt_2[,jj])

      var_2 <- var_2 + ( epsilon_1 %*% t(epsilon_2) )

      count <- count + 1
    }
    var_2 <- var_2 / count

    traza_2 <- dot_integrate(r = x$fparam, v = diag(var_2))

    sup_corr <- sup_cov / ( sqrt(traza_1)*sqrt(traza_2) )

    # Check - L2 surface norm
    vector_PACF[lag_PACF] <-
      dot_integrate(
        dot_integrate_col(v=t(sup_corr^2),r=x$fparam),
        r=x$fparam)

    # vector_PACF[lag_PACF] <-
    #   sqrt( .obtain_suface_L2_norm(x$fparam, list(Lag0 = sup_corr)) )
  }

  plt <- .plot_FACF(rho = vector_PACF, SWN = FACF$SWN_bound, WWN=NULL, ylab='PACF', ...)
  if(figure){
    plt
  }

  invisible( list('pacfs' = vector_PACF, 'SWN' = FACF$SWN_bound,
                  'WWN_bound'=NULL, 'plot'=plt) )
}


#' Estimate distribution of the fACF under the iid hypothesis using MC method
#'
#' Estimate the distribution of the autocorrelation function under the
#'  hypothesis of strong functional white noise. This function uses a
#'  Monte Carlo method to estimate the distribution.
#'
#' @inheritParams acf.dfts
#' @param autocovSurface An \eqn{(m x m)} matrix with the discretized
#' values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained
#' by calling the function \code{autocovariance()}.
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
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item Y <- generate_brownian_bridge(N=200, v=30) \\
#'          autocovSurface <- autocovariance(Y,0:lag.max) \\
#'          matindex <- .obtain_suface_L2_norm(v,autocovSurface)[-1] \\
#'          MC_dist <- .estimate_iid_distr_MC(Y,autocovSurface,matindex)
#'  }
#'
#' @keywords internal
#' @noRd
.estimate_iid_distr_MC <-
  function(x, autocovSurface, matindex, nsims= 10000){
    x <- dfts(x)

    # # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
    # # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
    # means <- matrix(rep(rowMeans(x$data),ncol(x$data)), nrow=nrow(x$data))
    # l <- .obtain_autocov_eigenvalues(x$data - means,x$fparam)
    l <- .obtain_autocov_eigenvalues(center(x))

    neig <- length(l)
    Reig <- rep(0,nsims)

    for(jj in 1:neig){
      for(kk in 1:neig){
        Reig <- Reig +
          stats::rchisq(n = nsims,df = 1)*l[jj]*l[kk]
      }
    }
    #Reig=Reig/nrow(Y)
    Reig <- Reig / ncol(x$data)

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
#' @inheritParams acf.dfts
#' @param epsilon Value used to determine how many eigenvalues will be returned.
#'   The eigenvalues \eqn{\lambda_{j}/\lambda_{1} > \code{epsilon}} will be
#'   returned. By default \code{epsilon = 0.0001}.
#'
#' @return A vector containing the \eqn{k} eigenvalues
#' greater than \code{epsilon}.
#'
#' @examples
#' #N <- 100
#' #v <- seq(from = 0, to = 1, length.out = 10)
#' #sig <- 2
#' #Y <- generate_brownian_bridge(N, v, sig)
#' #lambda <- .obtain_autocov_eigenvalues(x = Y)
#'
#' @keywords internal
#' @noRd
.obtain_autocov_eigenvalues <- function(x, epsilon = 0.0001){
  x <- dfts(x)

  nobs <- ncol(x$data) # nt
  res <- nrow(x$data) # nv

  # w(i,j) = integral of the product of the curves i and j
  W <- matrix(0,nrow = nobs,ncol = nobs)

  for(ii in 1:nobs){
    mat.aux <- matrix(rep(x$data[,ii],each = nobs),
                      nrow = nobs, ncol = res)*t(x$data)
    for(jj in 1:nobs){
      W[ii,jj] <- dot_integrate(r = x$fparam, v = mat.aux[jj,])
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
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item Y <- generate_brownian_bridge(N=200, v=30) \\
#'          autocovSurface <- autocovariance(Y,0:lag.max) \\
#'          matindex <- .obtain_suface_L2_norm(v,autocovSurface)[-1] \\
#'          MC_dist <- .estimate_iid_distr_Imhof(Y,autocovSurface,matindex)
#'  }
#'
#' @keywords internal
#' @noRd
.estimate_iid_distr_Imhof <- function(x, autocovs, l2norms){
  x <- dfts(x)

  if(!requireNamespace('CompQuadForm')) stop('Install `CompQuadForm` to run Imhof.')

  # # mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
  # # l <- obtain_autocov_eigenvalues(v,Y - mat.means)
  # means <- matrix(rep(rowMeans(x$data),ncol(x$data)), nrow=nrow(x$data))
  # l <- .obtain_autocov_eigenvalues(x$data - means,x$fparam)
  l <- .obtain_autocov_eigenvalues(center(x))

  nl <- length(l)
  Reig <- .estimate_iid_distr_MC(x, autocovs, l2norms)$Reig
  ex <- seq(from = 0, to = max(Reig), length.out = 250)
  x <- ex * ncol(x$data)
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
#'  distribution obtained by calling the function \code{acf.dfts}.
#' @param WWN The upper prediction bound for the weak white noise iid
#'  distribution obtained by calling the function \code{acf.dfts}.
#' @param ... Further arguments passed to the \code{plot} function.
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item fACF <- acf.dfts(x = generate_brownian_bridge(100), lag.max = 15,
#'                           alpha=0.95, figure = FALSE) \\
#'          .plot_FACF(rho = fACF$acfs,SWN = fACF$SWN_bound, WWN = fACF$WWN_bound)
#'  }
#'
#' @keywords internal
#' @noRd
.plot_FACF <- function(rho, SWN, WWN, ...){
  # Save preferences in case of an error
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  # Define suitable lwd for plotting
  lag.max <- length(rho)

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  arguments1 <- list()
  if(!"xlab" %in% names(arguments)){
    arguments1$xlab <- "Lag"
  }else{
    arguments1$xlab <- arguments$xlab
  }
  if(!"ylab" %in% names(arguments)){
    arguments1$ylab <- "ACF"
  }else{
    arguments1$ylab <- arguments$ylab
  }
  if(!"ylim" %in% names(arguments)){
    arguments1$ylim <- c(0, min(max(rho)*1.5,1))
  }else{
    arguments1$ylim <- arguments$ylim
  }
  if(!"lwd"  %in% names(arguments)){
    arguments1$lwd <- 2#lwd_1
  }else{
    arguments1$lwd <- arguments$lwd
  }
  if(!"las"  %in% names(arguments)){
    arguments1$las <- 1
  }else{
    arguments1$las <- arguments$las
  }
  if(!"lend" %in% names(arguments)){
    arguments1$lend <- 2
  }else{
    arguments1$lend <- arguments$lend
  }
  if(!"yaxs" %in% names(arguments)){
    arguments1$yaxs <- "i"
  }else{
    arguments1$yaxs <- arguments$yaxs
  }
  if(!"xaxs" %in% names(arguments)){
    arguments1$xaxs <- "i"
  }else{
    arguments1$xaxs <- arguments$xaxs
  }
  if(!"main" %in% names(arguments)){
    arguments1$main <- ""
  }else{
    arguments1$main <- arguments$main
  }
  if(!"xlim" %in% names(arguments)){
    arguments1$xlim <- c(0, length(rho)+1)
  }else{
    arguments1$xlim <- arguments$xlim
  }
  if(!"cex.axis" %in% names(arguments)){
    arguments1$cex.axis <- 1.5
  }else{
    arguments1$cex.axis <- arguments$cex.axis
  }
  if(!"cex.lab" %in% names(arguments)){
    arguments1$cex.lab <- 2.5
  }else{
    arguments1$cex.lab <- arguments$cex.lab
  }
  if(!"mar" %in% names(arguments)) graphics::par(mar=c(5,6,4,1)+.1)

  arguments1$x <- seq(1, length(rho), by = 1)
  arguments1$y <- rho
  arguments1$type <- "h"
  # arguments$xlim <- c(0,1)
  #arguments$ylim <- c(0,1.0)

  do.call(graphics::plot, arguments1)
  #grid(lty = 1)
  do.call(graphics::lines, arguments1)
  graphics::lines(x = arguments1$x,
                  y = arguments1$y,
                  type = arguments1$type,
                  col = 'black',#"lightgrey",
                  # lwd = 2,
                  # lwd = arguments1$lwd - 2,
                  lend = 2)
  blue_col <- "#0073C2FF"
  graphics::abline(h = SWN, col = blue_col, lty = 2)# lwd = 4, lty = 2)
  graphics::lines(x = arguments1$x, y = WWN, col = 'red', lty = 2)# lwd = 4, lty = 2)
  # graphics::legend(
  #   x = "topleft",
  #   legend = c(paste("i.i.d. bound (",ci*100," % conf.)",sep="")),
  #   col = blue_col,
  #   lty = 2,
  #   lwd = 4)
  graphics::box()
}

#' Fit an ARH(p) to a given functional time series
#'
#' Fit an \eqn{ARH(p)} model to a given functional time series. The fitted
#'   model is based on the model proposed in (Aue et al, 2015), first
#'   decomposing the original functional observations into a vector time series
#'   of \code{n_pcs} FPCA scores, and then fitting a vector autoregressive
#'   model of order \eqn{p} (\eqn{VAR(p)}) to the time series of the scores.
#'   Once fitted, the Karhunen-Loeve expansion is used to re-transform the
#'   fitted values into functional observations.
#'
#' @inheritParams pacf.dfts
#' @param p Numeric value specifying the order
#' of the functional autoregressive
#' model to be fitted.
#' @param show_varprop Logical. If \code{show_varprop = TRUE}, a plot of
#'   the proportion of variance explained by the first \code{n_pcs} functional
#'   principal components will be shown. By default \code{show_varprop = TRUE}.
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item mod <- .fit_ARHp_FPCA(x = generate_brownian_motion, p = 1, n_pcs = 5)
#'  }
#'
#' @references Aue, A., Norinho, D. D., Hormann, S. (2015).
#'  \emph{On the Prediction of Stationary Functional Time Series}
#'  Journal of the American Statistical Association,
#'  110, 378--392.
#'
#' @keywords internal
#' @noRd
.fit_ARHp_FPCA <- function(x, p, n_pcs, show_varprop = TRUE){
  x <- dfts(x)

  nobs <- ncol(x$data) #dt <- nrow(y)
  res <- nrow(x$data) #dv <- length(v)

  # Step 1: FPCA decomposition of the curves
  pca <- stats::princomp(t(x$data))
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
       fitted_vals = fitted_vals, x = x)
}
