#' PCA Characteristic Function Change Point Analysis
#'
#' This implements a principle component based characteristic change method.
#'
#' @inheritParams change
#'
#' @return List with pvalue, location, statistic, and simulations
#'
#' @references Huskova, M., & Meintanis, S.G. (2006). Change Point Analysis
#'  based on Empirical Characteristic Functions. Metrika, 63, 145-168.
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item res <- .change_pca_distribution(generate_brownian_motion(20,
#'                                          v=seq(0,1,length.out=10)), TVE=0.2)
#'  }
.change_pca_distribution <- function(X, statistic='Tn', critical='resample',
                                     TVE = 0.95, gam = 0.5, M = 200,
                                     blocksize = 1, resample_blocks = 'separate',
                                     replace = TRUE) {

  tmp <- .pca_distribution_statistic(X, statistic, TVE, gam, TRUE)
  stat <- tmp[1]
  location <- tmp[2]

  if(critical=='resample'){
    simulations <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                              type = resample_blocks, replace = replace,
                              fn = .pca_distribution_statistic,
                              statistic=statistic, TVE = TVE, gam=gam)
  }

  list('pvalue' = sum(stat <= simulations) / M,
       'location' = location)#,
       # 'statistic' = stat,
       # 'simulations'=as.vector(simulations))
}

#' Compute Test Statistic for PCA Distribution Change
#'
#' @inheritParams .change_pca_distribution
#' @param location Boolean for if location should be returned
#'
#' @returns Test statistic and potentially the location of the change
#'
#' @noRd
#' @keywords internal
.pca_distribution_statistic <- function(X, statistic, TVE, gam, location=FALSE){

  X_pca <- pca(dfts(X), TVE=TVE)
  n <- ncol(X)

  kappa <- matrix(nrow=length(X_pca$sdev),ncol=nrow(X_pca$x)-1)

  # Params for loop
  t_vals <- seq(0, 1, length.out=length(X_pca$x[,1]) )
  ws <- .w(t_vals)
  ks <- 1:(n-1)
  weights <- ((ks * (n - ks)) / n^2)^gam * ((ks * (n - ks)) / n)

  for(i in 1:length(X_pca$sdev)){
    Y <- X_pca$x[,i]

    dat <- exp(complex(real = 0, imaginary = 1) * t_vals %*% t(Y) )

    ## This commented code is a slow version of below
    # internals <- sapply(ks, function(k, dat, ws){
    #   abs(rowMeans(dat[,1:k,drop=FALSE]) - rowMeans(dat[,(k+1):n,drop=FALSE]) )^2 * ws
    # },dat=dat,ws=ws)

    cmean <- t(apply(dat,MARGIN = 1,cumsum)) /
      matrix(1:ncol(dat),nrow=length(t_vals),ncol=n,byrow = TRUE)
    cmean1 <- t(apply(dat[,n:2,drop=FALSE],MARGIN = 1,cumsum))[,(n-1):1] /
      matrix((n-1):1,nrow=length(t_vals),ncol=n-1,byrow = TRUE)
    internals <- abs(cmean[,-n] - cmean1)^2 * ws

    kappa[i,] <- weights * dot_integrate_col(internals)
  }
  if(length(X_pca$sdev)>1){
    Qnk <-  diag( 1/n * ( t(kappa) %*% diag(1/X_pca$sdev^2) %*% kappa ) )
  }else{
    Qnk <-  diag( 1/n * ( t(kappa) %*% (1/X_pca$sdev^2) %*% kappa ) )
  }

  if(statistic=='Tn'){
    stat <- dot_integrate(Qnk)
  }else if(statistic=='Mn'){
    stat <- max(Qnk)
  }

  if(!location) return(stat)

  location <- which.max(Qnk)

  c(stat, location)
}

# #' Compute Test Statistic for PCA Distribution Change
# #'
# #' @inheritParams .change_pca_distribution
# #' @param location Boolean for if location should be returned
# #'
# #' @returns Test statistic and potentially the location of the change
# #'
# #' @noRd
# #' @keywords internal
# .pca_distribution_statistic <- function(X, statistic, TVE, gam, location=FALSE){
#
#   X_pca <- pca(dfts(X), TVE=TVE)
#   n <- ncol(X)
#
#   kappa <- matrix(nrow=length(X_pca$sdev),ncol=nrow(X_pca$x)-1)
#   for(i in 1:length(X_pca$sdev)){
#     for (k in 1:(n - 1)) {
#       kappa[i,k] <- ((k * (n - k)) / n^2)^gam * ((k * (n - k)) / n) *
#         stats::integrate(
#           function(t, Y, k) {
#             abs(.phi_k(Y, t, k) - .phi_k0(Y, t, k))^2 * .w(t)
#           },
#           lower = 0, upper = 1, Y = X_pca$x[,i], k = k
#         )[[1]]
#     }
#   }
#   if(length(X_pca$sdev)>1){
#     Qnk <-  diag( 1/n * ( t(kappa) %*% diag(1/X_pca$sdev^2) %*% kappa ) )
#   }else{
#     Qnk <-  diag( 1/n * ( t(kappa) %*% (1/X_pca$sdev^2) %*% kappa ) )
#   }
#
#   if(statistic=='Tn'){
#     stat <- dot_integrate(Qnk)
#   }else if(statistic=='Mn'){
#     stat <- max(Qnk)
#   }
#
#   if(!location) return(stat)
#
#   location <- which.max(Qnk)
#
#   c(stat, location)
# }

# #' Compute Phi_k for scalar detection
# #'
# #' @inheritParams .change_pca_distribution
# #' @param t Time dimension
# #' @param k potential change
# #'
# #' @return Numeric
# #'
# #' @noRd
# #' @keywords internal
# .phi_k <- function(Y, t, k) {
#   sumVal <- 0
#   for (j in 1:k) {
#     sumVal <- sumVal +
#       exp(complex(real = 0, imaginary = 1) * t * Y[j])
#   }
#
#   1 / k * sumVal
# }
#
#
# #' Compute Phi_k0 for scalar detection
# #'
# #' @inheritParams .phi_k
# #'
# #' @return Numeric
# #'
# #' @noRd
# #' @keywords internal
# .phi_k0 <- function(Y, t, k) {
#   n <- length(Y)
#   sumVal <- 0
#
#   for (j in (k + 1):n) {
#     sumVal <- sumVal +
#       exp(complex(real = 0, imaginary = 1) * t * Y[j])
#   }
#
#   1 / (n - k) * sumVal
# }

#' Weighting Scheme for scalar change point detection
#'
#' @inheritParams .phi_k
#' @param a (Optional) Weight to consider. Default is 1.
#'
#' @return Numeric weight
#'
#' @noRd
#' @keywords internal
.w <- function(t, a = 1) {
  exp(-a * abs(t))
}
