

#' Estimate the Covariance Matrix
#'
#' This (internal) function computes the covariance matrix of some FD data.
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param x (Optional) Vector of locations of observations. Default is equally
#'     spaced observations on (0,1) with the same number of observations as FDs
#' @param M (Optional) Integer indicating the number of vectors used to create
#'     each value in the covariance matrix. Default is 25.
#' @param h (Optional) Integer indicating amount of lag to consider. Default is
#'     3.
#' @param K (Optional) Function for the kernel function to use. Default is the
#'     barlett_kernel.
#' @param W (Optional) Data.frame of vectors for measuring space. If NULL then
#'     then a Guassian measure is generated. Default is NULL.
#'
#' @return Data.frame for covariance based on Gaussian measure and given data X
#'
#' @examples
#' # This is an internal function, see usage in computeMethod
.estimCovMat <- function(X, x=seq(0,1,length.out=nrow(X)),
                         M=25, h=3, K=bartlett_kernel, W=NULL){
  # Setup random vectors
  if(is.null(W)){
    W <- computeSpaceMeasuringVectors(M,"BM",X)
  }
  if(ncol(W)!=M) stop('Number of Vectors not M')

  D11 <- D12 <- D21 <- D22 <- data.frame(matrix(ncol=M,nrow=M))
  funs <- c(cos,sin)
  for(i in 1:M){
    for(j in i:M){
      D11[i,j] <- D11[j,i] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[1][[1]],v=W[,i],
                lfunp=funs[1][[1]],vp=W[,j])
      D12[i,j] <- D21[j,i] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[1][[1]],v=W[,i],lfunp=funs[2][[1]],vp=W[,j])
      if(j>i){
        D12[j,i] <- D21[i,j] <- .estimD(K=K, h=h, X=X,
                                        lfun=funs[1][[1]], v=W[,j],
                                        lfunp=funs[2][[1]], vp=W[,i])
      }

      D22[i,j] <- D22[j,i] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[2][[1]],v=W[,i],lfunp=funs[2][[1]],vp=W[,j])
    }
  }

  tmp1 <- x %*% t(x)
  minDF <- matrix(1,nrow=length(x),ncol=length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      minDF[i,j] <- (min(x[i], x[j]) - tmp1[i,j])
    }
  }

  # Mults value from minDF to D11, then next value of minDF to D11, ...
  # rbind(cbind(fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D11)),
  #             fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D12))),
  #       cbind(fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D21)),
  #             fastmatrix::kronecker.prod(as.matrix(minDF),as.matrix(D22))))
  fastmatrix::kronecker.prod(
    as.matrix(minDF),
    rbind(cbind(as.matrix(D11),as.matrix(D12)),
          cbind(as.matrix(D21),as.matrix(D22))))

}


#' Estimate Long-run Covariance (D) matrix
#'
#' This (internal) function
#'
#' @param K Function for the kernel function to use.
#' @param h Integer indicating amount of lag to consider.
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#' @param lfunp Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param vp Numeric vector from L2 space with same number of observations as
#'     FD observation.
#'
#' @return Data.frame with autocovariance value based on data given
#'
#' @examples
#' # This is an internal function, see usage in .estimCovMat
.estimD <- function(K,h,X,lfun,v,lfunp,vp){

  iters <- (1-ncol(X)):(ncol(X)-1)

  # Move so pass in function values to avoid recomputation later
  fVals <- as.numeric(.estimf(X,lfun,v))
  fpVals <- as.numeric(.estimf(X,lfunp,vp))

  mean1 <- mean(fVals)
  mean2 <- mean(fpVals)

  values <- sapply(iters, function(k,K,h,X1,fVals,fpVals,mean1,mean2){
    Kval <- K(k,h)

    ifelse(Kval==0,
           0,
           Kval * .estimGamma(k=k,  X=X1,
                              fVals=fVals, fpVals=fpVals,
                              mean1=mean1, mean2=mean2))
  },K=K, h=h, X1=X, fVals=fVals, fpVals=fpVals,
  mean1=mean(fVals),mean2=mean(fpVals))

  sum(values) / ncol(X)
}


#' Estimate Autocovariance (gamma) Function
#'
#' This (internal) function computes the autocovariance (gamma) function.
#'
#' @param k Integer indicating FD object to consider
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#' @param lfunp Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param vp Numeric vector from L2 space with same number of observations as
#'     FD observation.
#'
#' @return Numeric autocovariance value for given data
#'
#' @examples
#' # This is an internal function, see usage in .estimD
.estimGamma <- function(k, X, fVals, fpVals, mean1, mean2){
  tmp <- 0

  if(k >=0){
    rs <- 1:(ncol(X)-k)
  }else{
    rs <- (1-k):ncol(X)
  }

  sum(.estimR(rs, fVals, mean1)*.estimR(rs+k, fpVals ,mean2))
}


#' Estimate R-hat
#'
#' This (internal) function estimates R-hat, that is lfun(<X_r,v>) - Y where Y
#'     demeans the data, typically mean(lfun(<X,v>)).
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param r Integer indicating the FD object of interest in the data X
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#' @param meanVal (Optional) Value indicating the meanVal. This can be is
#'     computed beforehand for speed. Default of NA computes mean(f-hat)
#'
#' @return Numeric value of R-hat for fd object of interest
#'
#' @examples
#' # This is an internal function, see usage in .estimGamma
.estimR <- function(r, fVals, meanVal=NA){
  if(is.na(meanVal)) meanVal <- mean(fVals)

  fVals[r] - meanVal
}


#' Estimate f-hat
#'
#' This (internal) function computes estimate of f-hat, that is function(<X,v>).
#'
#' @param Xr Numeric vector/data.frame indicating FD observation(s).
#' @param lfun Function, typically cos or sin. This is important when considering
#'     real or imaginary part.
#' @param v Numeric vector from L2 space with same number of observations as
#'     FD observation.
#'
#' @return Numeric value(s) indicating function value for each FD object given
#'
#' @examples
#' # This is an internal function, see usage in .estimR
.estimf <- function(Xr,lfun,v){
  lfun(t(Xr) %*% v)
}
