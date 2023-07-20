
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
#'     spaced observations on (0, 1) with the same number of observations as FDs
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
#' \dontrun{
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
#' }
detect_changepoint <- function(X, nSims=100, x=seq(0,1,length.out=ncol(X)),
                               M=25, h=3, K=bartlett_kernel, silent=FALSE){
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
      complex(imaginary = 1) * mvnorms[1,MJ+1:MJ]

    # Estimate value
    # gamProcess[i] <- sum(abs(gamVals)^2)/MJ
    gamProcess[i] <- .approx_int(abs(gamVals)^2)/nrow(X)
  }

  list('pval'=1-stats::ecdf(gamProcess)(value),
       'gamProcess'=gamProcess,
       'value'=value)
}

#' Estimate null and detect change point based on a single covar matrix
#'
#' @param X XXXXXX
#' @param nSims XXXXXX
#' @param x XXXXXX
#' @param h XXXXXX
#' @param K XXXXXX
#' @param space XXXXXX
#' @param silent XXXXXX
#' @param TN_M XXXXXX
#' @param Cov_M XXXXXX
#'
#' @return XXXXXX
#' @export
#'
#' @examples
#' \dontrun{
#' X <- generate_data_fd(ns = c(100,100),
#'                       eigsList = list(c(3,2,1,0.5),c(3,2,1,0.5)),
#'                       basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                                        fda::create.bspline.basis(nbasis=4, norder=4)),
#'                       meansList = c(0,1),
#'                       distsArray = c('Normal'),
#'                       evals = seq(0,1,0.05),
#'                       kappasArray = c(0.5),silent = TRUE)
#' cp_res <- detect_changepoint_singleCov(X, nSims=500, x=seq(0,1,length.out=20),
#'                                        h=3, K=bartlett_kernel, silent=FALSE)
#'
#' X1 <- generate_data_fd(ns = c(200),
#'                        eigsList = list(c(3,2,1,0.5)),
#'                        basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'                        meansList = c(0),
#'                        distsArray = c('Normal'),
#'                        evals = seq(0,1,0.05),
#'                        kappasArray = c(0),silent = TRUE)
#' nocp_res <- detect_changepoint_singleCov(X1, nSims=500, x=seq(0,1,length.out=20),
#'                                          h=0, K=bartlett_kernel, silent=FALSE)
#' X3 <- generate_data_fd(ns = c(50,250),
#'                        eigsList = list(c(3,2,1,0.5),c(30,1)),
#'                        basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                                         fda::create.bspline.basis(nbasis=2, norder=2)),
#'                        meansList = c(0,0),
#'                        distsArray = c('Normal'), kappasArray = c(0.5),
#'                        evals = seq(0,1,0.05), silent = TRUE)
#' cp_res3 <- detect_changepoint_singleCov(X, nSims=500,
#'                                         x=seq(0,1,length.out=20))
#' }
detect_changepoint_singleCov <- function(X, nSims=2000, x=seq(0,1,length.out=40),
                                         h=3, K=bartlett_kernel, space='BM',
                                         silent=FALSE, TN_M=10000, Cov_M=75){
  # Source Cpp File to speed up matrix computations
  # Rcpp::sourceCpp("R/matrixMult.cpp")

  # Determine Number of Iterations
  val_Tn <- compute_Tn(X,M=TN_M, space=space)

  # Generate Noise
  set.seed(123)
  W <- computeSpaceMeasuringVectors(Cov_M,space,X)

  # Variables
  MJ <- Cov_M * length(x)

  # Compute Gamma Matrix
  covMat <- .estimCovMat(X,x,Cov_M,h,K,W)
  covMat <- round(covMat,15)
  covMat_svd <- La.svd(covMat)
  sqrtD <- diag(sqrt(covMat_svd$d))
  sqrtD[is.na(sqrtD)] <- 0
  # sqrtMat <- eigenMapMatMult(
  #   eigenMapMatMult(covMat_svd$u, sqrtD),t(covMat_svd$v))
  sqrtMat <- Rfast::mat.mult(
    Rfast::mat.mult(covMat_svd$u,sqrtD),covMat_svd$vt)
  rm(W,covMat,covMat_svd,sqrtD)

  gamProcess <- c()
  nIters <- nSims/100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, lx){
    # (After trans + mult) Rows are iid MNV
    mvnorms <- sapply(1:100,function(m,x){stats::rnorm(x)},x=2*MJ)
    mvnorms <- Rfast::mat.mult(sqrtMat, mvnorms)

    gamVals <- mvnorms[1:MJ,] + complex(imaginary = 1)*mvnorms[MJ+1:MJ,]

    # Estimate value
    apply(abs(t(gamVals))^2,MARGIN = 1, .approx_int) / nrow(X)#nrow(X)
  }, MJ=MJ, sqrtMat=sqrtMat, lx=length(x))

  gamProcess <- as.vector(gamProcess)

  list('pval'=1-stats::ecdf(gamProcess)(val_Tn),
       'gamProcess'=gamProcess,
       'value'=val_Tn)
}


#' Estimate the Covariance Matrix
#'
#' This (internal) function computes the covariance matrix of some FD data.
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param x (Optional) Vector of locations of observations. Default is equally
#'     spaced observations on (0, 1) with the same number of observations as FDs
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
#' @noRd
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

  D11 <- D12_tD21 <- D22 <- data.frame(matrix(ncol=M,nrow=M))
  funs <- c(cos,sin)
  for(i in 1:M){
    for(j in i:M){
      D11[i,j] <- D11[j,i] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[1][[1]],v=W[,i],
                lfunp=funs[1][[1]],vp=W[,j])
      D12_tD21[i,j] <-
        .estimD(K=K,h=h,X=X,
                lfun=funs[1][[1]],v=W[,i],lfunp=funs[2][[1]],vp=W[,j])
      if(j>i){
        D12_tD21[j,i] <- .estimD(K=K, h=h, X=X,
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
    rbind(cbind(as.matrix(D11),as.matrix(D12_tD21)),
          cbind(as.matrix(t(D12_tD21)),as.matrix(D22))))

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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
#'
#' @examples
#' # This is an internal function, see usage in .estimR
.estimf <- function(Xr,lfun,v){
  lfun(t(Xr) %*% v)
}
