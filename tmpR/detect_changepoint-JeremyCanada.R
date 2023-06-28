
#' Detect Change Point
#'
#' This function detects change points in functional data. This is done through
#'     simulating the null distribution
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param nSims (Optional) Integer indicating the number of realizations of
#'     Gaussian processes to compute. Default is 100
#' @param x (Optional) Vector of locations of observations. Default is equally
#'     spaced observations on (0,1) with the same number of observations as FDs
#' @param M (Optional) Integer indicating the number of vectors used to create
#'     each value in the covariance matrix. Default is 25.
#' @param h (Optional) Integer indicating amount of lag to consider. Default is
#'     3.
#' @param K (Optional) Function for the kernel function to use. Default is the
#'     barlett_kernel.
#' @param silent (Optional) Boolean indicating it the output should be supressed.
#'     Default is FALSE.
#'
#' @return List
#'     1. pval: pvalue based on the data and estimated Gaussian processes
#'     2. gamProcess: Vector of estimated test statistics based on data
#'     3. value: Test statistic for data
#' @export
#'
#' @examples
#' X <- generate_data_fd(ns = c(100),
#'            eigsList = list(c(3,2,1,0.5)),
#'            basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'            meansList = c(0),
#'            distsArray = c('Normal'),
#'            evals = seq(0,1,0.05),
#'            kappasArray = c(0))
#' tmp <- detect_changepoint(X)
#' tmp$pval
#'
#' X1 <- generate_data_fd(ns = c(50,50),
#'            eigsList = list(c(3,2,1,0.5),
#'                            c(1.5,1,0.5,0.25)),
#'            basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                             fda::create.bspline.basis(nbasis=4, norder=4)),
#'            meansList = c(0,0),
#'            distsArray = c('Normal','Binomial'),
#'            evals = seq(0,1,0.05),
#'            kappasArray = c(0,0))
#' tmp1 <- detect_changepoint(X1)
#' tmp1$pval
detect_changepoint <- function(X, nSims=100, x=seq(0,1,length.out=ncol(X)),
                               M=25, h=3, K=bartlett_kernel, silent=F){
  # Setup Vars
  gamProcess <- rep(NA,nSims)
  value <- compute_Tn(X,M=M)
  MJ <- M * length(x)

  for(i in 1:nSims){
    if(!silent) cat(paste0(i,', '))
    # Compute Gamma Matrix
    covMat <- .estimCovMat(X,x,M,h,K)

    #mvnorms <- mvtnorm::rmvnorm(1,sigma=covMat,method='svd')
    mvnorms <- suppressWarnings(mvtnorm::rmvnorm(1,sigma=covMat))

    gamVals <- mvnorms[1,1:MJ] +
      complex(imaginary = 1)*mvnorms[1,MJ+1:MJ]

    # Estimate value
    gamProcess[i] <-
      sum(abs(gamVals[-c(MJ-0:(length(x)-1))])^2)/MJ
  }

  list('pval'=1-ecdf(gamProcess)(value), 'gamProcess'=gamProcess, 'value'=value)
}

#' Estimate null and detect change point based on a single covar matrix
#'
#' @param X
#' @param nSims
#' @param x
#' @param h
#' @param K
#' @param silent
#' @param maxM
#' @param ratio
#'
#' @return
#' @export
#'
#' @examples
#' X <- generate_data_fd(ns = c(100,100),
#'                       eigsList = list(c(3,2,1,0.5),c(3,2,1,0.5)),
#'                       basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                                        fda::create.bspline.basis(nbasis=4, norder=4)),
#'                       meansList = c(0,1),
#'                       distsArray = c('Normal'),
#'                       evals = seq(0,1,0.05),
#'                       kappasArray = c(0.5),silent = T)
#' cp_res <- detect_changepoint_singleCov(X, nSims=500, x=seq(0,1,length.out=20),
#'                                        h=3, K=bartlett_kernel, silent=F)
#'
#' X1 <- generate_data_fd(ns = c(200),
#'                        eigsList = list(c(3,2,1,0.5)),
#'                        basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'                        meansList = c(0),
#'                        distsArray = c('Normal'),
#'                        evals = seq(0,1,0.05),
#'                        kappasArray = c(0),silent = T)
#' nocp_res <- detect_changepoint_singleCov(X1, nSims=500, x=seq(0,1,length.out=20),
#'                                          h=1, K=bartlett_kernel, silent=F)
#' X3 <- generate_data_fd(ns = c(50,250),
#'                        eigsList = list(c(3,2,1,0.5),c(30,1)),
#'                        basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                                         fda::create.bspline.basis(nbasis=2, norder=2)),
#'                        meansList = c(0,0),
#'                        distsArray = c('Normal'), kappasArray = c(0.5),
#'                        evals = seq(0,1,0.05), silent = T)
#' cp_res3 <- detect_changepoint_singleCov(X, nSims=500,
#'                                         x=seq(0,1,length.out=20),
#'                                         h=3, K=bartlett_kernel, silent=F,
#'                                         maxM=250, ratio=0.05)
detect_changepoint_singleCov <- function(X, nSims=500, x=seq(0,1,length.out=20),
                                         h=3, K=bartlett_kernel, space='BM',
                                         silent=F,maxM=250, ratio=0.05){
  # Determine Number of Iterations
  done <- FALSE
  M <- 90
  prevValues <- c(compute_Tn(X,M=M), compute_Tn(X,M=M))

  while(!done && M<maxM){
    M <- M + 10
    values <- c(compute_Tn(X,M=M),compute_Tn(X,M=M))

    if(abs(1-values[1]/prevValues[1])<=ratio &&
       abs(1-values[2]/prevValues[2])<=ratio &&
       abs(1-values[1]/values[2])<=ratio)
      done <- TRUE
    prevValues <- values
  }

  # Generate Noise
  W <- computeSpaceMeasuringVectors(M,space,X)
  #data.frame(matrix(nrow=nrow(X),ncol=M))
  #for(i in 1:M){
  #  W[,i] <- sde::BM(N = nrow(X)-1)
  #}

  # Variables
  MJ <- M * length(x)

  # Compute Gamma Matrix
  covMat <- .estimCovMat(X,x,M,h,K,W)
  covMat_svd <- svd(covMat) # covMat_svd$u %*% diag(covMat_svd$d) %*% t(covMat_svd$v)
  sqrtD <- sqrt(diag(covMat_svd$d))
  sqrtD[is.na(sqrtD)] <- 0
  sqrtMat <- covMat_svd$u %*% sqrtD %*% t(covMat_svd$v)
  rm(W,covMat,covMat_svd,sqrtD)

  gamProcess <- c()
  nIters <- nSims/100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ,sqrtMat,lx){
    # (After trans + mult) Rows are iid MNV
    mvnorms <- t(sapply(1:100,function(m,x){rnorm(x)},x=2*MJ))
    mvnorms <- mvnorms %*% sqrtMat

    gamVals <- mvnorms[,1:MJ] + complex(imaginary = 1)*mvnorms[,MJ+1:MJ]

    # Estimate value
    rowSums(abs(gamVals[,-c(MJ-0:(lx-1))])^2)/MJ

  }, MJ=MJ, sqrtMat=sqrtMat, lx=length(x))

  gamProcess <- as.vector(gamProcess)
  # for(i in 1:nIters){
  #   # (After trans + mult) Rows are iid MNV
  #   mvnorms <- t(sapply(1:100,function(m,x){rnorm(x)},x=2*MJ))
  #   mvnorms <- mvnorms %*% sqrtMat
  #
  #   gamVals <- mvnorms[,1:MJ] + complex(imaginary = 1)*mvnorms[,MJ+1:MJ]
  #
  #   # Estimate value
  #   gamProcess <- c(gamProcess,
  #                   rowSums(abs(gamVals[,-c(MJ-0:(length(x)-1))])^2)/MJ)
  # }

  list('pval'=1-ecdf(gamProcess)(values[1]), 'pval2'=1-ecdf(gamProcess)(values[1]),
       'gamProcess'=gamProcess, 'value'=values[1], 'value2'=values[2],
       'M'=M)
}


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

  if(is.null(W)){
    W <- data.frame(matrix(nrow=nrow(X),ncol=M))
    for(i in 1:M){
      W[,i] <- sde::BM(N = nrow(X)-1)
    }
  }
  if(ncol(W)!=M) stop('Number of Vectors not M')

  D11 <- D12 <- D21 <- D22 <- data.frame(matrix(ncol=M,nrow=M))
  funs <- c(cos,sin)
  for(i in 1:M){
    for(j in i:M){
      D11[i,j] <- D11[j,i] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[1][[1]],v=W[,i],lfunp=funs[1][[1]],vp=W[,j])
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
  tmp2 <- minDF <- matrix(x,nrow=length(x),ncol=length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      minDF[i,j] <- (min(tmp2[i,j], t(tmp2)[i,j]) - tmp1[i,j])
    }
  }

  rbind(cbind(kronecker(as.matrix(minDF),as.matrix(D11)),
              kronecker(as.matrix(minDF),as.matrix(D12))),
        cbind(kronecker(as.matrix(minDF),as.matrix(D21)),
              kronecker(as.matrix(minDF),as.matrix(D22))))
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

  sumVal <- 0
  KVals <- sapply((1-ncol(X)):(ncol(X)-1), K, h=h)

  for(k in (1-ncol(X)):(ncol(X)-1)){
    if(KVals[k+ncol(X)]!=0){
      sumVal <- sumVal + KVals[k+ncol(X)] *
        .estimGamma(k=k,X=X,lfun=lfun,v=v,lfunp=lfunp,vp=vp)
    }
  }

  sumVal
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
.estimGamma <- function(k,X,lfun,v,lfunp,vp){
  tmp <- 0

  if(k >=0){
    mean1 <- mean(.estimf(X,lfun,v))
    mean2 <- mean(.estimf(X,lfunp,vp))
    for(r in 1:(ncol(X)-k)){
      tmp <- tmp + .estimR(X,r,lfun,v,mean1) * .estimR(X,r+k,lfunp,vp,mean2)
    }
  }else{
    mean1 <- mean(.estimf(X,lfun,v))
    mean2 <- mean(.estimf(X,lfunp,vp))
    for(r in (1-k):ncol(X)){
      tmp <- tmp + .estimR(X,r,lfun,v,mean1) * .estimR(X,r+k,lfunp,vp,mean2)
    }
  }

  1/ncol(X) * tmp
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
.estimR <- function(X,r,lfun,v, meanVal=NA){
  if(is.na(meanVal))
    meanVal <- mean(.estimf(X,lfun,v))

  .estimf(X[,r],lfun,v) - meanVal
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
