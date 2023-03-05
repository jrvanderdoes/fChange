#' Compute Spacing Measuring Vectors
#'
#' This function is used to compute vectors from spaces to measure the function
#'     data objects.
#'
#' @param M Integer for the number of vectors to generate
#' @param space String for the space of interest. Options are 'BM', 'OU' and 'PC'
#' @param X Data of interest. Used to know evaluation amount of vectors
#'
#' @return Data.frame with columns of vectors describing the space
#' @export
#'
#' @examples
#' X <- generate_data_fd(ns = c(50),
#'                       eigsList = list(c(3,2,1,0.5)),
#'                       basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'                       meansList = 0,
#'                       distsArray = c('Normal'),
#'                       evals = seq(0,1,0.05),
#'                       kappasArray = c(0))
#'
#' computeSpaceMeasuringVectors(10,'BM', X)
#' computeSpaceMeasuringVectors(10,'OU', X)
#' computeSpaceMeasuringVectors(10,'PC', X)
computeSpaceMeasuringVectors <- function(M, space, X){
  if(space=='BM'){
    W <- as.data.frame(sapply(rep(0,M),sde::BM, N=nrow(X)-1))
  } else if(space=='PC'){
    pComps <- prcomp(X, center=F, scale=F)
    W <- as.data.frame(sapply(rep(nrow(X),M),
                              function(x,pcs){
                                rowSums(rnorm(ncol(pcs))*pcs)
                              },pcs=pComps$x))
  } else if(space=='OU'){
    x <- seq(0,1, length.out=nrow(X))
    covMat <- (matrix(1,ncol=length(x),nrow=length(x)))
    covMat[,1] <- covMat[1,] <- exp(abs(x-x[1]))
    for(i in 2:(length(x)-1)){

      covMat[,i] <- covMat[i,] <- c(covMat[i,c(1:(i-1))],
                                    exp(abs(x-x[i]))[-c(1:(i-1))])
    }
    W <- as.data.frame(sapply(rep(nrow(X),M),
                              function(x,covMat){
                                covMat %*% rnorm(x)
                              },covMat=covMat))

    # W <- as.data.frame(sapply(rep(nrow(X),M),sde::rsOU,theta=c(0,0,1)))
  } else{
    stop('Error: Sorry only BM, PC, or OU processes allowed')
  }

  W
}
