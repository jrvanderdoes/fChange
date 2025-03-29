#' Testing changes in the trace of the covariance operator in functional data
#'
#' This function tests and detects changes in the trace of the covariance operator.
#'
#' @details This function dates and detects changes in trace of the covariance
#'  function. This can be interpreted as the changes in the total variation
#'  of the the functional data. Trace is defined as the infinite sum of
#'  the eigenvalues of the covariance operator and for the sake of
#'  implementation purpose, the sum is truncated up to the total number
#'  of basis functions that defines the functional data at hand. The critical
#'  values are approximated via \code{M} Monte Carlo simulations.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param changes If not NULL then the data is centered considering the changes.
#' @param M Number of monte carlo simulations used to get the critical values. The default value is \code{M=1000}
#'  value is \code{h=2}.
#'
#' @return
#' \item{\code{pvalue}}{
#'  Approximate p value for testing whether there is a significant change in the desired eigenvalue of the covariance operator
#' }
#' \item{\code{change}}{
#'  Estimated change location
#' }
#' \item{\code{trace_before}}{
#'  Estimated trace before the change
#' }
#' \item{\code{trace_after}}{
#'  Estimated trace after the change
#' }
#'
#' @noRd
#' @keywords internal
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2020). Structural break
#'  analysis for spectrum and trace of covariance operators. Environmetrics
#'  (London, Ont.), 31(1)
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .change_trace(generate_brownian_motion(200))
#'  }
.change_trace <- function(X, changes = NULL, M = 1000,
                          statistic = 'Tn', critical='simulation',
                          blocksize=1, replace=TRUE, type='separate'){
  X <- center(dfts(X), changes=changes)

  N <- ncol(X$data)
  D <- nrow(X$data)

  tmp <- .trace_statistic(X$data, statistic, location=TRUE)
  Sn <- tmp[1]
  k_star <- tmp[2]

  if(critical=='simulation'){
    # values_sim <- sapply(1:M, function(k,N) abs(sde::BBridge(0, 0, 0, 1, N-1)), N=N )
    values_sim <- abs(generate_brownian_bridge(M,seq(0,1,length.out=N))$data)

    if(statistic=='Tn'){
      Values <- dot_integrate_col(values_sim)
    }else if(statistic=='Mn'){
      Values <- apply(values_sim, MARGIN = 2, max)
    }
  } else if(critical=='resample'){
    Values <- .bootstrap(X = X$data, blocksize = blocksize, M = M,
                         type = type, replace = replace, fn = .trace_statistic,
                         statistic=statistic)
  }
  # Tn <- rep(0,N)
  # for (k in 1:N) {
  #   T_x <- sum(Xi[1:k])/N
  #   Tn[k] <- (1/sigma) * abs(T_x - (k)/N * T_1)
  # }
  # Sn <- max(Tn)
  p <- sum(Sn <= Values) / M # Compute p-value
  # k_star <- min(which.max(Tn))

  # data1_pca <- pca(dfts(X$data[,1:k_star]),TVE=1)
  # data2_pca <- pca(dfts(X$data[,(k_star+1):N]),TVE=1)
  # tr_before <- sum(data1_pca$sdev)
  # tr_after <- sum(data2_pca$sdev)
  list('pvalue' = p,
       'location' = k_star)#,
       # 'statistic' = Sn,
       # 'simulations' = Values
       # )
}


#' Trace Test Statistic
#'
#' function for test statistic in trace change point detection
#'
#' @inheritParams .change_trace
#' @param location Boolean if location of change should be returned
#'
#' @returns Test statistic and (optional) location of change
#'
#' @noRd
#' @keywords internal
.trace_statistic <- function(X, statistic, location=FALSE){
  N <- ncol(X)

  Cov_op <- .partial_cov(X, 1)
  lambda <- Cov_op$eigen_val
  T_1 <- sum(lambda)
  Xi <- sapply(1:N, function(i) X[,i] %*% X[,i])
  sigma_sq <- sandwich::lrvar(Xi, prewhite = FALSE)
  sigma <- sqrt(sigma_sq)

  Tn <- rep(0,N)
  for (k in 1:N) {
    T_x <- sum(Xi[1:k])/N
    Tn[k] <- (1/sigma) * abs(T_x - k/N * T_1)
  }

  if(statistic=='Tn'){
    Sn <- dot_integrate(Tn)
  }else if(statistic=='Mn'){
    Sn <- max(Tn)
  }

  if(!location) return(Sn)

  c(Sn, min(which.max(Tn)) )
}


#' Partial Sample Estimates of the Covariance Function in Functional Data Analysis
#'
#' This function computes the partial sum estimate of the covariance function
#'  in functional data analysis. It also generate the eigenvalues and
#'  eigenfunctions of the parial sum estimate, along with the coefficient
#'  matrix of the estimated covariance function.
#'
#' @details This function simply estimates the covariance function based on the partial sum of the centered functional observations.
#'  The length of the sum is determined by \code{x}, and when \code{x=1} this estimate corresponds to the regular covariance function
#'  estimates using the whole sample.
#'
#' @inheritParams .change_trace
#' @param x Fraction of the sample size that the partial sum is computed.
#'  This input must be in (0,1]. The default is \code{x=1} which corresponds
#'  to the regular covariance function estimate in functional data analysis.
#'
#' @keywords internal
#' @noRd
#'
#' @return
#'  \item{\code{eigen_val}}{
#'    Eigenvalues of the partial sum estimate of the covariance function
#'  }
#'  \item{\code{eigen_fun}}{
#'    Eigenfunctions of the partial sum estimate of the covariance function
#'  }
#'  \item{\code{coef_matrix}}{
#'    Coefficient matrix of the partial sum estimate of the covariance function
#'  }
#'
#' @seealso \code{\link{pca.fd}}
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .partial_cov(electricity)
#'  }
.partial_cov <- function(X, x = NULL){
  X <- dfts(X)
  if (is.null(x)) x <- 1
  if (x>1 | x<=0) stop("x should be in (0,1]")

  n <- ncol(X$data)
  D <- nrow(X$data)
  k <- floor(n*x)

  X <- center(X)
  C <- matrix(0, D, D)
  for (i in 1:k){
    C <- C + X$data[,i] %*% t(X$data[,i])
  }
  E <- eigen(C/n)
  e_val <- E$values
  e_fun <- E$vectors # fd(E$vectors, fdobj$basis)

  list(eigen_val = E$values,
       eigen_fun = E$vectors,
       coef_matrix = C/n)
}

