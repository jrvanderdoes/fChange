#' Fully Functional Mean Change Point Analysis
#'
#' This function tests whether there is a significant change in the mean
#'     function of the functional data, and it will give an estimate for the
#'     location of the change. The procedure is based on the standard L-2 norm
#'     and hence does not depend on any dimension reduction technique such as
#'     fPCA.
#'
#' @param data A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param statistic String \code{Tn} or \code{Mn}, for integrated or maximum test statistic
#' @param M (Optional) Number of Monte Carlo simulations used to get the critical
#'     values. The default value is \code{M=1000}
#' @param  h (Optional) The window parameter parameter for the estimation of the
#'     long run covariance kernel. The default value is \code{h=0}, i.e., it
#'     assumes iid data
#' @param K (Optional) Function indicating the Kernel to use if h>0
#'
#' @noRd
#' @keywords internal
#'
#' @return If inc.pval is false, Numeric of CP location or NA, otherwise
#'     location or NA and the p-value
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2018). Detecting and dating structural
#'  breaks in functional data without dimension reduction. Journal of the Royal
#'  Statistical Society. Series B, Statistical Methodology, 80(3), 509-529.
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .change_mean(generate_brownian_motion(500), M = 250)
#'    \item .change_mean(electricity, M = 250)
#'  }
.change_mean <- function(data, statistic=c('Mn','Tn'), critical='simulation',
                         M = 1000, h = 0, K = bartlett_kernel,
                         blocksize=1, type = 'separate',replace = TRUE) {
  statistic <- .verify_input(statistic, c('Tn','Mn'))
  data <- dfts(data)
  data <- center(data)

  ## Estimate eigenvalues (lambda_i, 1<=i<=d)
  Ceps <- long_run_covariance(data, h, K) # data1
  lambda <- eigen(Ceps)$values

  tmp <- .mean_statistic(data$data, statistic,location=TRUE)
  stat <- tmp[1]
  k.star <- tmp[2]

  if(critical=='simulation'){
    values_sim <- sapply(1:M, function(k, lambda, n, statistic)
      .asymp_dist(n, lambda, statistic),
      lambda = lambda, n = ncol(data), statistic = statistic
    )
  } else if(critical=='resample'){
    values_sim <- .bootstrap(X = data, blocksize = blocksize, M = M,
                             type = type, replace = replace, fn = .mean_statistic,
                             statistic=statistic)
  }

  p <- sum(stat <= values_sim) / M


  list('pvalue' = p,
       'location' = k.star)#,
       # 'change_fun' = mean(data$data[,1:k.star])-mean(data$data[,(k.star+1):n]),
       # 'statistic' = stat, 'simulations' = values_sim)
}


#' Statistic for Mean Change
#'
#' @inheritParams .change_mean
#' @param location Boolean if location should also be returned
#'
#' @returns Either (1) the test statistic value or (2) the test statistic value
#'  and the estimated location
#'
#' @noRd
#' @keywords internal
.mean_statistic <- function(data, statistic, location=FALSE){
  n <- ncol(data)

  Sn2 <- rep(0, n)
  # CUSUM
  for (k in 1:n) {
    # TODO:: add weights
    # Sn2[k] <- compute_mean_stat(data,k=k,weight=0.5)
    Sn2[k] <- sum((rowSums(data[, 1:k,drop=FALSE]) -
                     (k / n) * rowSums(data))^2)
  }
  Sn2 <- Sn2 / n

  # Compute Test statistic
  if(statistic=='Tn'){
    stat <- dot_integrate(Sn2)
  }else if(statistic=='Mn'){
    stat <- max(Sn2, na.rm = TRUE)
  }else{
    stop('Set `statistic` to `Tn` or `Mn`.',call. = FALSE)
  }

  if(!location) return(stat)

  c(stat, min(which(Sn2 == max(Sn2,na.rm = TRUE))))
}


#' Asymptotic Distribution
#'
#' Define and simulate asymptotic distribution based on Brownian bridge and
#'     eigenvalues.
#'
#' @param n Integer indicating the number of observations for the Brownian bridge
#' @param lambda Vector of numerics indicating eigenvalues
#'
#' @return Numeric indicating max value of the Brownian bridge
#'
#' @noRd
#' @keywords internal
.asymp_dist <- function(n, lambda, statistic='Mn') {
  BridgeLam <- matrix(0, length(lambda), n)

  BridgeLam <-
    t(generate_brownian_bridge(length(lambda),v=seq(0,1,length.out=n))$data^2 %*%
    diag(lambda))

  if(statistic=='Tn'){
    threshold <- dot_integrate(colSums(BridgeLam))

  } else if(statistic=='Mn'){
    threshold <- max(colSums(BridgeLam))
  }

  threshold
}


# #' Compute CUSUM Statistic for Mean Change
# #'
# #' This function is used to compute the CUSUM statistic for a mean change of
# #'     FD observations.
# #'
# #' @param data dfts object or numeric data.frame with evaled points on rows and fd objects in columns
# #' @param k Numeric indicating the change point of interest
# #' @param ... Unused, just for use in other functions
# #'
# #' @return Numeric indicating the CUSUM value at the given K value
# #'
# #' @noRd
# #' @keywords internal
# .compute_mean_stat <- function(data, k, weight = 0.5) {
#   data <- dfts(data)
#   n <- ncol(data$data)
#   stop('Not tested with weights')
#   normalizer <- ( (k/n) * ((n-k) / n) )^(-weight)
#   sum(normalizer * (rowSums(data$data[, 1:k,drop=FALSE]) -
#          (k / n) * rowSums(data$data))^2) / n
# }

