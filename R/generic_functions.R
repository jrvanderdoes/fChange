#' Generic Centering of Data
#'
#' @param object Object for computation of centering
#' @param CPs Change points for centering individual sections.
#' @param type String of \code{mean} or \code{median} for method of centering.
#' @param ... Unused
#'
#' @name center
#'
#' @return Centered data of the same format as the given data
#' @export
#'
#' @examples
#' center(1:10)
NULL

#' Generic pca of Data
#'
#'  Build on \code{prcomp}.
#'
#' @param object Object for computation of principle components analysis
#' @param ... Additional parameters to function based on data
#'
#' @name pca
#'
#' @return PCA data
#' @export
#'
#' @examples
#' pca(1:10)
NULL

#' Generic sd/var of Data
#'
#' @param object Object for computation of standard deviation or variance of
#'  a given data set.
#' @param ... Additional parameters to function based on data
#'
#' @name sd
#'
#' @return Numeric(s) to explain the standard deviation / variance
#' @export
#'
#' @examples
#' sd(1:10)
#' var(1:10)
NULL


#############################################


#' @rdname center
#'
#' @seealso
#'  [center.default()], [center.data.frame()],
#'  [center.matrix()], [center.funts()]
#'
#' @export
center <- function(object, CPs=NULL, type='mean', ...) UseMethod("center")
#' @rdname center
#' @export
center.default <- function(object, CPs=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    method <- mean
  }else if(type=='median'){
    method <- stats::median
  }


  if(is.null(CPs)){
    object <- object - method(object)
  }else{
    CPs_use <- unique(c(0, CPs, length(object)))
    CPs_use <- CPs_use[order(CPs_use)]
    for(i in 1:(length(CPs_use)-1)){
      object[(CPs_use[i]+1):CPs_use[i+1]] <-
        object[(CPs_use[i]+1):CPs_use[i+1]] -
        method(object[(CPs_use[i]+1):CPs_use[i+1]])
    }
  }

  object
}
#' @rdname center
#' @export
center.data.frame <- function(object, CPs=NULL, type='mean', ...) { object <- as.matrix(object); NextMethod("center") }
#' @rdname center
#' @export
center.matrix <- function(object, CPs=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    row_method <- rowMeans
  }else if(type=='median'){
    row_method <- function(y){apply(y, MARGIN = 1, stats::median)}
  }


  if(is.null(CPs)){
    object <- object - row_method(object)
  }else{
    CPs_use <- unique(c(0, CPs, ncol(object)))
    CPs_use <- CPs_use[order(CPs_use)]
    for(i in 1:(length(CPs_use)-1)){
      object[,(CPs_use[i]+1):CPs_use[i+1]] <-
        object[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] -
        row_method(object[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE])
    }
  }

  object
}
#' @rdname center
#' @export
center.funts <- function(object, CPs=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    row_method <- rowMeans
  }else if(type=='median'){
    row_method <- function(y){apply(y, MARGIN = 1, stats::median)}
  }


  if(is.null(CPs)){
    object$data <- object$data - row_method(object$data)
  }else{
    CPs_use <- unique(c(0, CPs, ncol(object$data)))
    CPs_use <- CPs_use[order(CPs_use)]
    for(i in 1:(length(CPs_use)-1)){
      object$data[,(CPs_use[i]+1):CPs_use[i+1]] <-
        object$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] -
        row_method(object$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE])
    }
  }

  object
}



#' Generic function for PCA
#'
#' @param object Data to compute principal component analysis on
#' @param TVE Numeric in \[0,1\] for the total variance explained. Can use this to
#'  select only the required components.
#' @param ... Additional information for \code{prcomp}.
#'
#' @return Principal component data
#'
#' @export
pca <- function(object, TVE = 1,...) UseMethod("pca")
#' @rdname pca
#' @export
pca.default <- function(object, ...) { stats::prcomp(object, ...) }


#' @rdname pca
#'
#' @export
#'
#' @examples
#' pca(funts(electricity))
pca.funts <- function(object, TVE = 1, ...){
  if(TVE > 1 || TVE < 0) stop('TVE must be in [0,1].',call. = FALSE)

  pc <- stats::prcomp(x=t(object$data), ...)
  if(TVE==1){
    min_pc <- length(pc$sdev)
  } else{
    min_pc <- min( which(cumsum(pc$sdev^2)/sum(pc$sdev^2) > TVE) )
  }

  ## TODO:: Check out
  new_rot <- pc$rotation[,1:min_pc,drop=FALSE]*sqrt(nrow(object$data))

  # scores <- matrix(NA, nrow=ncol(x$data),ncol = min_pc)
  # for(i in 1:ncol(x$data)){
  #   for(l in 1:min_pc){
  #     scores[i,l] <- dot_integrate(x$data[,i] %*% t(new_rot[,l]))
  #   }
  # }
  scores <- (t(object$data-pc$center) %*% new_rot)[,1:min_pc,drop=FALSE]/nrow(object$data)

  # res = fda::pca.fd(fda::Data2fd(X$data))
  # pc$sdev^2
  # res$values[1:2]
  # View(pc$x); View(res$scores)
  # View(pc$rotation); View(fda::eval.fd(seq(0,1,length.out=30),res$harmonics))

  list(sdev = pc$sdev[1:min_pc]/sqrt(nrow(object$data)),
       rotation = new_rot,
       center = pc$center,
       scale = pc$scale,
       x = scores,
       orig_pc = pc)
}


#' Generic Variance Standard Deviation
#'
#' @param object Data to compute on
#' @param ... Additional parameters based on the data
#'
#' @return Standard Deviation/Variance
#' @export
sd <- function(object, ...) UseMethod("sd")
#' @rdname sd
#' @export
sd.default <- function(object, ...) stats::sd(object, ...)

#' @rdname sd
#' @param type Character to specify if an operator ('op') or pointwise ('pw')
#'  calculation is desired
#'
#' @export
#'
#' @examples
#' sd(funts(electricity),type='pointwise')
sd.funts <- function(object, type='pointwise', ...) {
  type <- c('pointwise')[min(pmatch(type,c('pointwise')))]

  if(type=='pointwise'){
    return( apply(object$data, MARGIN = 1, sd) )
  } else{
    stop('Only type="pointwise" is implemented', call. = FALSE)
  }

}


#' @rdname sd
#' @export
var <- function(object, ...) UseMethod("var")
#' @rdname sd
#' @export
var.default <- function(object, ...) stats::var(object, ...)


#' @rdname sd
#'
#' @export
#'
#' @examples
#' var(funts(electricity),type='pointwise')
#' var(funts(electricity),type='operator')
var.funts <- function(object, type=c('operator','pointwise'), ...) {
  type <- c('operator','pointwise')[min(pmatch(type,c('operator','pointwise')))]

  if(type=='operator'){
    # autocov_approx_h(object$data,0)
    autocovariance(object$data,0)
  } else if(type=='pointwise'){
    apply(object$data, MARGIN = 1, var)
  } else{
    stop('Type must be "operator" or "pointwise"', call. = FALSE)
  }
}


#' Cumsum for funts object
#'
#' @param x funts object
#'
#' @return funts object with data as cumsum
#' @export
#'
#' @examples
#' cumsum(funts(electricity))
cumsum.funts <- function(x){
    funts(t(apply(x$data,MARGIN = 1,cumsum)),
          labels = x$labels,
          intraobs = x$intraobs)
}
