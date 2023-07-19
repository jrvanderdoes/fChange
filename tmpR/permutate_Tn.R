# Write and update documentation

#' Generalized Bootstrap and Permutation Methods Including Blocking
#'
#' Use this method for generalized resampling of a test statistic in
#'     determination of the proper cutoff value
#'
#' @param X A data frame of the evaluated fd data
#'     (rows are evaled points and columns are fd)
#' @param blockSize Numeric value indicating blocking sizes
#'
#' Value of 1 should be used for iid data while larger values can be used to
#'     account for dependence in the data
#'
#' @param fn Function that returns the test statistic value
#' @param iters Numeric value indicating number iterations used in bootstrapping
#' @param replace Boolean value indicating bootstrapping (T) or permutation (F)
#' @param alpha Numeric value in (0,1) to get the correct percentile
#' @param silent Boolean to silence output helpful to a user
#' @param ... Optional arguments passed into \code{fn}
#'
#' @return Numeric value indicating the cutoff value
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.fourier.basis(nbasis=2)),
#'     meansList = c(0,0.5),
#'     distsArray = c('Normal','Binomial'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0, 0.5))
#' # Permutation Method for Tn (this will get 1-alpha quantile of iters)
#' result <- permutating_Tn(X=data_KL,
#'     blockSize=ncol(data_KL)^(1/3),
#'     iters=1000, replace=FALSE)
permutate_Tn <- function(X, blockSize, iters = 1000,
                         W=computeSpaceMeasuringVectors(M=1000,space='BM',X=X),
                         replace=FALSE, alpha=0.05){

  .compute_Tn_CF <- function(CF){
    .IntVal <- function(CF){
      n <- nrow(CF)
      nx <- 0:n
      # Get nFhat - first value will be 0
      nfx <- t(data.frame(0,sapply(nx[-1],
                                function(x, CF){
                                  if(x==1){
                                    retVal <- CF[1,]
                                  } else{
                                    retVal <- colSums(CF[1:x,])
                                  }
                                  retVal}, CF=CF)))

      # Compute abs(Zn^2), Zn= sqrt(n)*(f(v,x)-floor(nx)/n*f(v,1))
      #   sqrt(n)/n = 1/sqrt(n)
      fVals <- apply(nfx, MARGIN = 2,
            function(val){
              t(abs(1/sqrt(n)*(val-nx/n*val[length(val)]))^2)
            }, simplify = TRUE)

      # LH Integration
      rowSums(fVals[,-ncol(fVals)]) * 1/n
    }

    1/ncol(CF) * sum(.IntVal(CF))
  }

  .compute_Loc <- function(CF){

    n <- nrow(CF)
    nx <- 0:n
    # Get nFhat - first value will be 0
    nfx <- t(data.frame(0,sapply(nx[-1],
                                 function(x, CF){
                                   if(x==1){
                                     retVal <- CF[1,]
                                   } else{
                                     retVal <- colSums(CF[1:x,])
                                   }
                                   retVal}, CF=CF)))

    # Compute abs(Zn^2), Zn= sqrt(n)*(f(v,x)-floor(nx)/n*f(v,1))
    #   sqrt(n)/n = 1/sqrt(n)
    fVals <- apply(nfx, MARGIN = 2,
                   function(val){
                     t(abs(1/sqrt(n)*(val-nx/n*val[length(val)]))^2)
                   }, simplify = TRUE)[-1,]

    list('loc'=which.max(1/nrow(fVals) *
                     rowSums(fVals)),
         'locPlot'=ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::geom_line(ggplot2::aes(x=nx[-1],
                                            y=1/nrow(fVals) * rowSums(fVals)))
         )
  }

  # Cols are cfs for one set (each ind)
  CF <- exp(complex(real=0,imaginary = 1) * (t(X) %*% as.matrix(W)))
  fullVal <- .compute_Tn_CF(CF)

  idxGroups <- .getChunks(1:nrow(CF), nrow(CF)/blockSize)
  idxs <- sapply(1:iters,function(i,m, indxs,replace){
    samps <- sample(1:m, replace = replace)
    unlist(indxs[samps], use.names = F)
  }, m=length(idxGroups),indxs=idxGroups, replace=replace)
  # If it is a matrix/dataframe this already even rows
  #     (will be with permutation test)
  if(!(is.matrix(idxs)| is.data.frame(idxs)))
    idxs <- .convertSamplesToDF(idxs)


  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter,X1){
                        .compute_Tn_CF(X1[stats::na.omit(loop_iter),]) },
                      X1=CF)

  append(list('cutoff'=stats::quantile(bssamples,
                                       probs = c(1-alpha))[[1]],
              'stat'=fullVal),
         .compute_Loc(CF))
}


#' Get Chunks
#'
#' This (internal) function splits the vector x into a chunksN number of
#'     subsegments. The values of x are kept in order (i.e. no scrambled).
#'
#' @param x Vector of values
#' @param chunksN Numeric indicating the number of chunks to split X into
#'
#' @return A list with chunksN items, each containing an similiar sized subset
#'     of the original vector
#'
#' @examples
#' .getChunks(1:100,1)
#' .getChunks(1:100,2)
#' .getChunks(1:100,5)
.getChunks <- function(x,chunksN) {
  if(chunksN<2)
    return(x)
  split(x, cut(x, chunksN, labels = FALSE))
}


#' Convert List of Samples into a Data Frame
#'
#' This (internal) function takes a list with differ length data.frames or
#'     vectors and pads them all to make a clean data.frame.
#'
#' @param data_list List of elements to be combined to a data.frame.
#'
#' @return Data.frame of the data in data_list
#'
#' @examples
#' # This is an internal function and will not be viewable to user. See
#' #     generalized_resampling.
.convertSamplesToDF <- function(data_list){

  m <- length(data_list)
  maxLen <- 0
  for(ii in 1:m){
    maxLen <- max(maxLen,length(data_list[[ii]]))
  }

  data_df <- data.frame(matrix(nrow=maxLen,ncol=m))

  for(ii in 1:length(data_list)){
    data_df[,ii] <- c(data_list[[ii]],
                      rep(NA,maxLen-length(data_list[[ii]])))
  }

  data_df
}
