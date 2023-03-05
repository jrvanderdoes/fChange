
#' Welsh Approximation to Tn
#'
#' This function approximates the Tn statistic using the Welsh approximation.
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param alpha (Optional) Numeric value in (0,1) for significance. Default is 0.05
#' @param TVal (Optional) Numeric value indicate the number of FDs. Default is the
#'     number of columns in X.
#' @param W (Optional) Data.frame for the functions to integrate against
#'     used for both muHat and sigmaHat2. Default of NULL will use Brownian
#'     motion.
#' @param W1 (Optional)  data.frame for the functions to integrate against
#'     used for sigmaHat2. Default of NULL will use Brownian motion.
#' @param M (Optional) Numeric value for number of functions to integrate against
#'     (ncol(W)). Default is 100.
#' @param h (Optional) Numeric value indicating the amount of blocking. Default
#'     is Tval^(1/3)
#' @param K (Optional) Function used for CHat, indicating Kernel. Default is the
#'     Bartlett Kernel.
#' @param ... Unused parameters, but added for use in other functions.
#'
#' @return Numeric value indicating cutoff from Welsh approximation
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,1),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' welsh_approximation(data_KL)
welsh_approximation <- function(X, alpha = 0.05, TVal = length(X[1,]),
                                W = NULL, W1 = NULL, M=100,
                                h = TVal^(1/3), K = bartlett_kernel, ...){
  if(is.null(W)){
    W <- as.data.frame(sapply(rep(0,M),sde::BM, N=length(X[,1])-1))
  }else{
    W <- as.data.frame(W)
  }
  if(is.null(W1)){
    W1 <- as.data.frame(sapply(rep(0,M),sde::BM, N=length(X[,1])-1))
  }else{
    W1 <- as.data.frame(W1)
  }

  muVals <- rep(0,M)
  sigma2Vals <- rep(0,M)

  for(i in 1:M){
    v <- W[,i]
    v1 <- W1[,i]

    CHat <- estimateLRV(X=X,v1=v,v2=v,TVal=TVal,h=h,K=K)
    CHat1 <- .estimateLRV(X=X,v1=v,v2=v1,TVal=TVal,h=h,K=K)

    muVals[i] <- CHat
    sigma2Vals[i] <- CHat1 * Conj(CHat1) # pairs
  }

  muHat <-  1/6 * mean(muVals)
  sigma2Hat <- (1/90) * 2 * mean(sigma2Vals)

  betaHat <- Re(sigma2Hat / (2*muHat))
  nuHat <- Re((2*muHat^2) / sigma2Hat)

  betaHat * qchisq(1-alpha, df=nuHat)
}


#' Estimate Long Run Variance
#'
#' This (internal) function estimates CHat
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param v1 Vector of length nrow(X) indicating a random vector in the space
#' @param v2 Vector of length nrow(X) indicating a random vector in the space
#' @param TVal Numeric value indicate the number of FDs.
#' @param h Numeric value indicating the amount of blocking.
#' @param K Function for the Kernel.
#'
#' @return Numeric indicated estimated CHat value
#'
#' @examples
#' # This is an internal function, see welsh_approximation for usage.
.estimateLRV <- function(X, v1, v2, TVal, h, K){

  eps1 <- exp(complex(real=0,imaginary = 1) * (t(X) %*% v1))
  #eps1Bar <- mean(eps1)#1/TVal * sum(eps1)

  eps2 <- exp(complex(real=0,imaginary = 1) * (t(X) %*% v2))
  #eps2Bar <- mean(eps2)#1/TVal * sum(eps2)

  gammaHat <- rep(NA, length((1-TVal):(TVal-1)))

  for(l in (1-TVal):(TVal-1)){
    gammaHat[TVal+l] <- .getAutocov(eps1=eps1, eps2=eps2,
                                  #eps1Bar=eps1Bar, eps2Bar=eps2Bar,
                                  l=l, TVal=TVal)
  }

  .computeLRVEstimate(h=h, TVal=TVal, gammaHat=gammaHat, K=K)
}


#' Estimate Autocovariance Estimates
#'
#' This (internal) function estimates autocovariance at single point
#'
#' @param eps1 Complex vector indicating exp(i, <X,v>)
#' @param eps2 Complex vector indicating exp(i, <X,v2>)
#' @param l Integer indicating lags
#' @param TVal Numeric value indicate the number of FDs.
#'
#' @return Numeric indicating autocovariance value at a particular lag
#'
#' @examples
#' # This is an internal function, see .estimateLRV for usage.
.getAutocov <- function(eps1, eps2, l, TVal){#eps1Bar, eps2Bar, l, TVal){
  gammaVal <- 0
  eps1Bar <- mean(eps1)
  eps2Bar <- mean(eps2)

  if(l >= 0){
    for(j in (1+l):TVal){
      gammaVal <- gammaVal +
        (eps1[j-l]-eps1Bar) * Conj(eps2[j]-eps2Bar)
    }
  } else{
    for(j in (1-l):TVal){
      gammaVal <- gammaVal +
        (eps1[j+l] - eps1Bar) * Conj(eps2[j]-eps2Bar) # Confirm this j+l swap
    }
  }

  1/TVal * gammaVal
}


#' Compute Long Run Variance
#'
#' This (internal) function computes long run variance estimates
#'
#' @param h Numeric value indicating the amount of blocking.
#' @param TVal Numeric value indicate the number of FDs.
#' @param gammaHat Vector of numerics for gamma estimates at each lag.
#' @param K Function for the Kernel.
#'
#' @return Numeric estimating LRV value
#'
#' @examples
#' # This is an internal function, see .estimateLRV for usage.
.computeLRVEstimate <- function(h, TVal, gammaHat, K){
  CVal <- 0

  for(l in (1-TVal):(TVal-1)){
    CVal <- CVal + K(l, h) * gammaHat[TVal+l]
  }

  CVal
}
