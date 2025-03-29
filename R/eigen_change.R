
#' Eigenvalue Changes - Detect eigenvalue changes in the Covariance Operator
#'
#' This function tests and detects changes individual or jointly in the
#'  eigenvalue(s) of the covariance operator.
#'
#' @details This function dates and detects changes in the joint eigenvalues
#'  that is defined by \code{d} of the covariance function. The critical
#'  values are approximated via \code{M} Monte Carlo simulations.
#'
#' @param X A functional data object of class '\code{fd}'
#' @param d Number of eigenvalues to include in testing.
#' @param h The window parameter for the estimation of the long run covariance matrix. The default
#'  value is \code{h=2}.
#' @param changes Vector to use in order to demean the data
#' @param test String if testing single or joint eigenvalues
#' @param M Number of monte carlo simulations used to get the critical values.
#' The default value is \code{M=1000}; however, simulations have shown this can
#' be greatly reduced. Consider this when using resample.
#' @param K Kernel function
#'
#' @noRd
#' @keywords internal
#'
#' @return
#' \item{\code{pvalue}}{
#'  Approximate p value for testing whether there is a significant change in the desired eigenvalue of the covariance operator
#' }
#' \item{\code{change}}{
#'  Estimated change location
#' }
#' \item{\code{eval_before}}{
#'  Estimated eigenvalues before the change
#' }
#' \item{\code{eval_after}}{
#'  Estimated eigenvalues after the change
#' }
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2020). Structural break
#'  analysis for spectrum and trace of covariance operators. Environmetrics
#'  (London, Ont.), 31(1)
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .change_eigen(electricity, 2, test='joint')
#'    \item .change_eigen(electricity, 2, test='individual')
#'
#'  }
.change_eigen <- function(X, d, h=2, changes = NULL,
                          statistic = c('Tn','Mn'),
                          test = c('joint', 'individual'),
                          critical = c('simulation','resample'),
                          M = 1000,
                          K=bartlett_kernel,
                          blocksize=1, type='separate', replace=TRUE){
  # Verification
  poss_tests <- c('joint', 'individual')
  test <- .verify_input(test,poss_tests)

  # Setup Data
  X <- center(dfts(X), changes=changes)

  N <- ncol(X$data)
  D <- nrow(X$data)

  # Eigen Change
  # Cov_op <- .partial_cov(X, 1)
  #
  # eig2 <- array(dim = c(D,D,D))
  # for(j in 1:D){
  #   eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
  # }
  # thetas <- matrix(nrow=N,ncol=D)
  # X2 <- array(dim = c(D,D,N))
  # for(i in 1:N){
  #   X2[,,i] <- X$data[,i] %*% t(X$data[,i]) - Cov_op$coef_matrix
  #
  #   for(j in 1:D){
  #     thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
  #   }
  # }

  if(test=='joint'){
    # Sigma_d <- long_run_covariance(X = t(thetas[,1:d]),h=h,K = K)

    # Values <- sapply(1:M,
    #                  function(i,N,d) .asymp_joint(N,d),
    #                  N=N, d=d)

    if(critical=='simulation'){

      if(statistic=='Tn'){
        Values <- sapply(1:M,
                         function(i,N,d) {
                           # tmp <- sapply(1:d,
                           #               function(dd,N) sde::BBridge(0, 0, 0, 1, N-1),
                           #               N=N)
                           tmp <- generate_brownian_bridge(d,v=seq(0,1,length.out=N))$data
                           dot_integrate(rowSums(tmp^2))
                         }, N=N, d=d)
      } else if(statistic=='Mn'){
        Values <- sapply(1:M,
                         function(i,N,d) {
                           # tmp <- sapply(1:d,
                           #               function(dd,N) sde::BBridge(0, 0, 0, 1, N-1),
                           #               N=N)
                           tmp <- generate_brownian_bridge(d,v=seq(0,1,length.out=N))$data
                           max(rowSums(tmp^2))
                         }, N=N, d=d)
      }

    } else if(critical=='resample'){
      Values <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                           type = type, replace = replace, fn = .eigen_joint_statistic,
                           statistic=statistic, d=d, h=h, K=K)
    }

    # Test statistic
    tmp <- .eigen_joint_statistic(X$data, statistic, d, h, K,  TRUE)
    Sn <- tmp[1]
    k_star <- tmp[2]
    p <- sum(Sn <= Values) / M # Compute p-value

    # before <- pca(X$data[,1:k_star])$sdev[1:d]
    # after <- pca(X$data[,(1+k_star):N])$sdev[1:d]

  } else if(test=='individual'){

    if(critical=='simulation'){

      if(statistic=='Tn'){
        # Values <- sapply(1:M, function(k)
        #   dot_integrate(sde::BBridge(0, 0, 0, 1, N-1)^2))
        Values <- dot_integrate_col(
          generate_brownian_bridge(M,v=seq(0,1,length.out=N))$data^2)
      } else if(statistic=='Mn'){
        # Values <- sapply(1:M, function(k)
        #   max(sde::BBridge(0, 0, 0, 1, N-1)^2))
        Values <- apply(
          generate_brownian_bridge(M,v=seq(0,1,length.out=N))$data^2,
          MARGIN=2, FUN = max)
      }

    } else if(critical=='resample'){
      Values <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                           type = type, replace = replace, fn = .single_eigen_statistic,
                           statistic=statistic, d=d, h=h, K=K)
    }

    # Test statistic
    tmp <- .single_eigen_statistic(X$data, statistic, d, h, K, TRUE)
    Sn <- tmp[1]
    k_star <- tmp[2]
    p <- sum(Sn <= Values) / M # Compute p-value


    # before <- pca(X$data[,1:k_star])$sdev[d]
    # after <- pca(X$data[,(1+k_star):N])$sdev[d]
  }else{
    stop('test must be "joint" or "individual".', call. = FALSE)
  }

  list('pvalue' = p,
       'location' = k_star)#,
       # 'statistic' = Sn,
       # 'simulations' = Values)
       # eval_before = before,
       # eval_after = after)
}



#' Eigen Change Test Statistic
#'
#' @inheritParams .change_eigen
#' @param X Matrix of data
#' @param location Boolean for if location should be returned as well
#'
#' @returns The test statistic and optional location
#'
#' @noRd
#' @keywords internal
.eigen_joint_statistic <- function(X, statistic, d, h, K, location=FALSE){

  N <- ncol(X)
  D <- nrow(X)

  Cov_op <- .partial_cov(X, 1)

  eig2 <- array(dim = c(D,D,D))
  for(j in 1:D){
    eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
  }
  thetas <- matrix(nrow=N,ncol=D)
  X2 <- array(dim = c(D,D,N))
  for(i in 1:N){
    X2[,,i] <- X[,i] %*% t(X[,i]) - Cov_op$coef_matrix

    for(j in 1:D){
      thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
    }
  }

  Sigma_d <- long_run_covariance(X = t(thetas[,1:d]), h = h, K = K)
  # Moore Penrose solve if non-invertable
  Sigma_d_inv <- tryCatch({
    solve(Sigma_d)
  }, error = function(e){
    eigs <- eigen( Sigma_d )
    eigs$vectors <- eigs$vectors * sqrt(D)
    eigs$values <- eigs$values / D

    K_eigs <- which.min(cumsum(eigs$values)/sum(eigs$values) < 0.95)

    tmp_inv <- matrix(0, nrow=nrow(Sigma_d), ncol=ncol(Sigma_d))
    for(i in 1:K_eigs){
      tmp_inv <- tmp_inv +
        ( 1 / eigs$values[i] ) * (eigs$vectors[,i] %o% eigs$vectors[,i] )
    }

    tmp_inv
  })

  Tns <-  rep(0,N)
  for(k in 1:N){
    lam <- .partial_cov(X, k/N)$eigen_val[1:d]
    kapa <- lam - (k/N) * Cov_op$eigen_val[1:d]
    Tns[k] <- N * t(kapa) %*% Sigma_d_inv %*% kapa
    # Tns[k] <- N * t(kapa) %*% solve(Sigma_d) %*% kapa
  }

  if(statistic=='Tn'){
    Sn <- dot_integrate(Tns)
  } else if(statistic=='Mn'){
    Sn <- max(Tns)
  }

  if(!location) return(Sn)

  c(Sn, min(which.max(Tns)))
}




#' Single Eigen Test Statistic
#'
#' @inheritParams .eigen_joint_statistic
#'
#' @returns Statistic or statistic and location for single eigen test
#'
#' @noRd
#' @keywords internal
.single_eigen_statistic <- function(X, statistic, d, h, K, location=FALSE){

  N <- ncol(X)
  D <- nrow(X)

  Cov_op <- .partial_cov(X, 1)

  eig2 <- array(dim = c(D,D,D))
  for(j in 1:D){
    eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
  }
  thetas <- matrix(nrow=N,ncol=D)
  X2 <- array(dim = c(D,D,N))
  for(i in 1:N){
    X2[,,i] <- X[,i] %*% t(X[,i]) - Cov_op$coef_matrix

    for(j in 1:D){
      thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
    }
  }

  Sigma_d <- long_run_covariance(X = t(thetas[,d]), h = h, K = K)

  Tns <- rep(0,N)
  for (k in 1:N){
    lam <- .partial_cov(X, k/N)$eigen_val[d]
    Tns[k] <- ( N*( lam - (k/N)*Cov_op$eigen_val[d] )^2 ) / Sigma_d
  }

  if(statistic=='Tn'){
    Sn <- dot_integrate(Tns)
  } else if(statistic=='Mn'){
    Sn <- max(Tns)
  }

  if(!location) return(Sn)

  c(Sn, min(which.max(Tns)))
}

# #' Asympotic Results for Joint
# #'
# #' Compute the test statistic threshold when joint
# #'
# #' @param N Length of data
# #' @param d Number of dimension
# #'
# #' @return Test statistic threshold
# #'
# #' @keywords internal
# #' @noRd
# .asymp_joint <- function(N, d){
#   WW <- matrix(0, d, N)
#   for(i in 1:d){
#     WW[i, ] <- sde::BBridge(0, 0, 0, 1, N-1)
#   }
#
#   max(colSums(WW^2))
# }
