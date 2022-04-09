
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
#' # Metric
#' compute_Tn(data_KL)
#' # Permutation Method for Tn
#' generalized_resampling(X=data_KL,
#'     blockSize=ncol(data_KL)^(1/3),
#'     fn=compute_Tn, iters=1000, replace=F)
generalized_resampling <- function(X, blockSize, fn, iters,
                                   replace=F, alpha=0.05, silent=F, ...){

  if(!silent){
    st<-Sys.time()
    fullVal <- fn(X,...)
    en<-Sys.time()

    cat(paste0('Estimated time: ',
               round(difftime(en,st,'units'='mins')[[1]]*iters,2),
               ' mins\n'))
  }

  n <- length(X[1,])
  idxGroups <- .getChunks(0:n, n/blockSize)
  idxs <- sapply(1:iters,function(i,m, indxs,replace){
    samps <- sample(1:m, replace = replace)
    unlist(indxs[samps])
  },m=length(idxGroups),indxs=idxGroups,replace=replace)
  # If it is a matrix/dataframe this already even rows
  #     (will be with permutation test)
  if(!(is.matrix(idxs)| is.data.frame(idxs)))
    idxs <- .convertSamplesToDF(idxs)


  bssamples <- sapply(as.data.frame(idxs),
                      function(loop_iter,fn,X1,...){ fn(X1[,na.omit(loop_iter)], ...) },
                      X1=X,fn=fn,...)

  quantile(bssamples,probs = c(1-alpha))[[1]]
}


.getChunks <- function(x,chunksN) split(x, cut(x, chunksN, labels = FALSE))

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
