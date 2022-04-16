
#' welsh Approximation to Tn
#'
#' @param X Numeric data.frame with evaled points on rows and fd objects in columns
#' @param alpha Numeric value in (0,1) for significance
#' @param TVal Numeric value indicate the number of FDs
#' @param W Optional data.frame for the functions to integrate against
#'     used for both muHat and sigmaHat2
#' @param W1 Optional data.frame for the functions to integrate against
#'     used for sigmaHat2
#' @param M Numeric value for number of functions to integrate against (ncol(W))
#' @param h Numeric value indicating the amount of blocking
#' @param K Function used for CHat
#'
#' @return Numeric value indicating cutoff
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
                                h = TVal^(1/3), K = bartlett_kernel){
  ## Code
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

    CHat <- .estimateCHat(X=X,v1=v,v2=v,TVal=TVal,h=h,K=K)
    CHat1 <- .estimateCHat(X=X,v1=v,v2=v1,TVal=TVal,h=h,K=K)

    muVals[i] <- CHat
    sigma2Vals[i] <- CHat1 * Conj(CHat1) # pairs
  }

  muHat <-  1/6 * (1/M * sum(muVals))
  sigma2Hat <- (1/90) * 2 * (1/M *  sum(sigma2Vals)) # Should M and T be the same??

  betaHat <- Re(sigma2Hat / (2*muHat))
  nuHat <- Re((2*muHat^2) / sigma2Hat)

  betaHat * qchisq(1-alpha, df=nuHat)
}

#' Bartlett Kernel
#'
#' One of the most obvious options for the K function in the Welsh-approximation
#'     is the Barlett Kernel.
#'
#' @param l A numeric value indicating current lag
#' @param h A numeric value indicating smallest lag with no influence
#'
#' @return A numeric value indicating value of Kernel at that point
#' @export
#'
#' @examples
#' bartlett_kernel(1,2)
#'
#' # Use in Welsh Approximation
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,1),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' welsh_approximation(data_KL,K=bartlett_kernel)
bartlett_kernel <- function(l, h){
  val <- 0
  if(abs(l) < h)
    val <- (1 - abs(l) / h)

  val
}


.estimateCHat <- function(X,v1,v2,TVal,h, K){
  ## Functions
  .getGamma <- function(eps1, eps2, eps1Bar, eps2Bar, l, TVal){
    gammaVal <- 0

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

  .getC <- function(h, TVal, gammaHat, K){
    CVal <- 0

    for(l in (1-TVal):(TVal-1)){
      CVal <- CVal + K(l, h) * gammaHat[TVal+l]
    }

    CVal
  }

  ## Code
  eps1 <- exp(complex(real=0,imaginary = 1) * (t(X) %*% v1))
  eps1Bar <- mean(eps1)#1/TVal * sum(eps1)

  eps2 <- exp(complex(real=0,imaginary = 1) * (t(X) %*% v2))
  eps2Bar <- mean(eps2)#1/TVal * sum(eps2)

  gammaHat <- rep(NA, length((1-TVal):(TVal-1)))

  for(l in (1-TVal):(TVal-1)){
    gammaHat[TVal+l] <- .getGamma(eps1=eps1, eps2=eps2,
                                  eps1Bar=eps1Bar, eps2Bar=eps2Bar,
                                  l=l, TVal=TVal)
  }
  CHat <- .getC(h=h, TVal=TVal, gammaHat=gammaHat, K=K)

  CHat
}


