
#' Change Point Analysis - Dimension Reduction
#'
#' This function tests whether there is a significant change in the mean
#'  function of functional data, and it gives an estimate of the location
#'  of the change. The procedure will reduce the dimension of the
#'  functional data using functional principal component analysis and will
#'  use the leading principal curves which explain \code{TVE} total variance
#'  to carry out the change point analysis.
#'
#' @details This functions performs structural break analysis for the functional
#'  data using an fPCA based initial dimension reduction. It is recommended
#'  that the dimension of the subspace, \code{d}, that the functional
#'  observations are projected onto should be selected based on TVE\.
#'
#' @inheritParams change
#' @param TVE Numeric for total variance explained to select the number of
#'  principle components
#' @param M Number of Monte Carlo simulations to get the critical values. The default value is \code{M=1000}
#' @param h The window parameter for the estimation of the long run covariance kernel. The default
#'  value is \code{h=0}, i.e., it assumes iid data
#' @param K Kernel function
#'
#' @return
#'\item{\code{pvalue}}{
#' An approximate p value for testing whether there is a significant change in the mean function
#'}
#'\item{\code{change}}{
#' Estimated change location
#'}
#'\item{\code{DataBefore}}{
#' Data before the estimated change
#'}
#'\item{\code{DataAfter}}{
#' Data after the estimated change
#'}
#'\item{\code{MeanBefore}}{
#' Mean function before the estimated change
#'}
#'\item{\code{MeanAfter}}{
#' Mean function after the estimated change
#'}
#'\item{\code{change_fun}}{
#' Estimated change function
#'}
#'
#' @references Berkes, I., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Detecting changes in the mean of functional observations}
#'  Journal of the Royal Statistical Society, Series B 71, 927-946
#' @references Aue, A., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Estimation of a change-point in the mean function of functional data}
#'  Journal of Multivariate Analysis 100, 2254-2269.
#'
#' @noRd
#' @keywords internal
#'
#' @examples
#' #res <- .change_pca_mean(generate_brownian_bridge(200,seq(0,1,length.out=10)))
#' #res1 <- .change_pca_mean(generate_ou(20,200))
#' #res2 <- .change_pca_mean(electricity)
.change_pca_mean <- function(X, statistic, critical, TVE=0.95,
                            M=1000, K=bartlett_kernel,
                            blocksize=1, resample_blocks='separate', replace = TRUE){
  X <- center(dfts(X))

  pca_X <- pca(dfts(X), TVE=TVE)
  d <- length(pca_X$sdev)

  D <- nrow(X$data)
  n <- ncol(X$data)

  # Test Statistic
  tmp <- .pca_mean_statistic(X, TVE, statistic, TRUE)
  stat <- tmp[1]
  location <- tmp[2]

  # Critical Values
  if(critical=='simulation'){

    if(statistic=='Tn'){

      simulations <- sapply(1:M, function(k,d,n){
        # B.Bridges <- rep(0, d)
        # for(j in 1:d){
        #   B.Bridges[j] <- dot_integrate( sde::BBridge(0,0,0,1,n-1)^2 )
        # }
        B.Bridges <-
          dot_integrate_col(generate_brownian_bridge(d,v=seq(0,1,length.out=n))$data^2)
        sum(B.Bridges)
      }, d=d,n=1000)
    } else if(statistic=='Mn'){

      simulations <- sapply(1:M, function(k,d,n){
        # B.Bridges <- matrix(NA, d,n)
        # for(j in 1:d){
        #   B.Bridges[j,] <-  sde::BBridge(0,0,0,1,n-1)^2
        # }
        B.Bridges <-
          t(generate_brownian_bridge(d,v=seq(0,1,length.out=n))$data^2)
        max(colSums(B.Bridges))
      }, d=d,n=D)
    }

  } else if(critical=='resample'){

    simulations <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                              type = resample_blocks, replace = replace,
                              fn = .pca_mean_statistic,
                              statistic=statistic, TVE = TVE)
  }

  list('pvalue' = sum(stat <= simulations) / M,
       'location' = location)#,
       # 'statistic' = stat,
       # 'simulations'=simulations)
}


#' Title
#'
#' @inheritParams .change_pca_mean
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param location Boolean if location should be returned or not
#'
#' @returns Test statistic and location (if requested)
#'
#' @noRd
#' @keywords internal
.pca_mean_statistic <- function(X, TVE, statistic, location=FALSE){

  # TODO:: Change to LRC
  pca_X <- pca(dfts(X), TVE=TVE)

  n <- ncol(X)
  d <- length(pca_X$sdev)

  eta.hat <- as.matrix(pca_X$x)

  ## Test Statistic
  kappa <-
    sapply(1:n,function(k, eta.hat, n){
      colSums(eta.hat[1:k, , drop=FALSE]) - k/n * colSums(eta.hat)
    },eta.hat=eta.hat, n=n)
  # If same names, only 1 pc and sapply flips it incorrectly
  if(length(unique(names(kappa)))==1) kappa <- t(kappa)

  if(length(pca_X$sdev)>1){
    Qnk <-  diag( 1/n * ( t(kappa) %*% diag(1/pca_X$sdev^2) %*% kappa ) )
  }else{
    Qnk <-  diag( 1/n * ( t(kappa) %*% (1/pca_X$sdev^2) %*% kappa ) )
  }

  if(statistic=='Tn'){
    # Snd_tmp <- rep(NA,d)
    # kappa1 <- matrix(ncol=n,nrow=d)
    # for(l in 1:d){
    #   inner_sums <- rep(0, n)
    #   for(k in 1:n){
    #     inner_sums[k] <- sum(eta.hat[1:k,l]) - k/n * sum(eta.hat[,l])
    #   }
    #   kappa1[l,] <- inner_sums
    #   Snd_tmp[l] <- 1/pca_X$sdev[l]^2 * sum(inner_sums^2)
    # }
    # stat <- 1/n^2 * sum(Snd_tmp)

    stat <- dot_integrate(Qnk)

  }else if(statistic=='Mn'){

    stat <- max(Qnk)

  }

  if(!location) return(stat)


  # kappa <-
  #   sapply(1:n,function(k, eta.hat, n){
  #     colSums(eta.hat[1:k, , drop=FALSE]) - k/n * colSums(eta.hat)
  #   },eta.hat=eta.hat, n=n)
  # Q_nk <- rep(NA,n)
  # for(k in 1:n){
  #   Q_nk[k] <- 1/n * ( t(kappa[,k]) %*% diag(1/pca_X$sdev^2) %*% kappa[,k] )
  # }
  location <- which.max(Qnk)

  c(stat, location)
}


# #' Compute S_n Test statistic for PCA change
# #'
# #' @param eta.hat XXX
# #' @param k XXX
# #' @param n XXX
# #'
# #' @return Numeric test statistic
# #'
# #' @keywords internal
# #' @noRd
# .S_n_pca <- function(eta.hat, k, n){
#   # normalizer = ((k/n) * ((n-k) / n))^(-0.5)
#   if(is.null(dim(eta.hat))){
#     eta.bar <- sum(eta.hat)/n
#     out <- sum(eta.hat[1:k]) - k*eta.bar
#   } else {
#     eta.bar <- as.matrix(rowSums(eta.hat)/n)
#     out <- rowSums(as.matrix(eta.hat[, 1:k])) - k*eta.bar
#   }
#
#   out #* normalizer
# }

# #' Asymptotic Threshold for mean PCA change
# #'
# #' @param N description
# #' @param d description
# #'
# #' @return Numeric asymptotic threshold
# #'
# #' @keywords internal
# #' @noRd
# .asymp_pca <- function(N, d){
#   B.Bridges <- matrix(0,nrow = d,ncol = N)
#   for(j in (1:d)){
#     B.Bridges[j,] <- sde::BBridge(0,0,0,1,N-1)^2
#   }
#   max(colSums(B.Bridges))
# }
