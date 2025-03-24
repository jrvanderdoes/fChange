#' Generic Centering of Data
#'
#' Center data by removing the mean or median. Defining changes allow for
#'  regional centering.
#'
#' @param object Object for computation of centering.
#' @param changes Change points for centering individual sections.
#' @param type String of \code{mean} or \code{median} for method of centering.
#' @param ... Parameters that may be fed into other versions of centering.
#'
#' @name center
#'
#' @return Centered data of the same format as the given data.
#' @export
#'
#' @examples
#' center(1:10)
NULL

#' Generic Function for Principal Component Analysis
#'
#' This is a generic function to call PCA on various objects. The default method
#'  uses [stats::prcomp()].
#'
#'
#' @param object Object for computation of principle components analysis.
#' @param TVE Numeric in \[0,1\] for the total variance explained, this determines
#'  the number of components and can be used for dimension reduction.
#' @param ... Additional parameters to extensions based on data. Often this is
#'  additional information for \code{prcomp}.
#'
#' @return Principal component data.
#' @export
#'
#' @name pca
#'
#'
#' @examples
#' pca(1:10)
NULL

#############################################


#' @rdname center
#'
#' @seealso
#'  [center.default()], [center.data.frame()],
#'  [center.matrix()], [center.dfts()]
#'
#' @export
center <- function(object, changes=NULL, type='mean', ...) UseMethod("center")
#' @rdname center
#' @export
center.default <- function(object, changes=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    method <- mean
  }else if(type=='median'){
    method <- stats::median
  }


  if(is.null(changes)){
    object <- object - method(object)
  }else{
    changes_use <- unique(c(0, changes, length(object)))
    changes_use <- changes_use[order(changes_use)]
    for(i in 1:(length(changes_use)-1)){
      object[(changes_use[i]+1):changes_use[i+1]] <-
        object[(changes_use[i]+1):changes_use[i+1]] -
        method(object[(changes_use[i]+1):changes_use[i+1]])
    }
  }

  object
}
#' @rdname center
#' @export
center.data.frame <- function(object, changes=NULL, type='mean', ...) { object <- as.matrix(object); NextMethod("center") }
#' @rdname center
#' @export
center.matrix <- function(object, changes=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    row_method <- rowMeans
  }else if(type=='median'){
    row_method <- function(y){apply(y, MARGIN = 1, stats::median)}
  }


  if(is.null(changes)){
    object <- object - row_method(object)
  }else{
    changes_use <- unique(c(0, changes, ncol(object)))
    changes_use <- changes_use[order(changes_use)]
    for(i in 1:(length(changes_use)-1)){
      object[,(changes_use[i]+1):changes_use[i+1]] <-
        object[,(changes_use[i]+1):changes_use[i+1],drop=FALSE] -
        row_method(object[,(changes_use[i]+1):changes_use[i+1],drop=FALSE])
    }
  }

  object
}
#' @rdname center
#' @export
center.dfts <- function(object, changes=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    row_method <- rowMeans
  }else if(type=='median'){
    row_method <- function(y){apply(y, MARGIN = 1, stats::median)}
  }


  if(is.null(changes)){
    object$data <- object$data - row_method(object$data)
  }else{
    changes_use <- unique(c(0, changes, ncol(object$data)))
    changes_use <- changes_use[order(changes_use)]
    for(i in 1:(length(changes_use)-1)){
      object$data[,(changes_use[i]+1):changes_use[i+1]] <-
        object$data[,(changes_use[i]+1):changes_use[i+1],drop=FALSE] -
        row_method(object$data[,(changes_use[i]+1):changes_use[i+1],drop=FALSE])
    }
  }

  object
}



#' @rdname pca
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
#' pca(electricity)
pca.dfts <- function(object, TVE = 1, ...){
  if(TVE > 1 || TVE < 0) stop('TVE must be in [0,1].',call. = FALSE)

  pc <- stats::prcomp(x=t(object$data), ...)
  if(TVE==1){
    min_pc <- length(pc$sdev)
  } else{
    min_pc <- min( which(cumsum(pc$sdev^2)/sum(pc$sdev^2) > TVE) )
  }

  ## TODO:: Check out
  new_rot <- pc$rotation[,1:min_pc,drop=FALSE] * sqrt(nrow(object$data))

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




#' Generic Function for Variance and Standard Deviation Computation
#'
#' Generic function to compute the variance and standard deviations. The default
#'  uses [stats::var()] and [stats::sd()].
#'
#' @param object Object for computation of standard deviation or variance of
#'  the given data set.
#' @param ... Additional parameters for the particular extensions.
#'
#' @name sdvar
#'
#' @return Numeric(s) to explain the standard deviation / variance
#' @export
#'
#' @seealso [stats::sd()], [stats::var()]
#'
#' @examples
#' sd(1:10)
#' var(1:10)
sd <- function(object, ...) UseMethod("sd")
#' @rdname sdvar
#' @export
sd.default <- function(object, ...) stats::sd(object, ...)

#' @rdname sdvar
#' @param type String to specify if an operator ('op') or pointwise ('pw')
#'  calculation is desired on the functional data.
#'
#' @export
#'
#' @examples
#' sd(electricity,type='pointwise')
sd.dfts <- function(object, type='pointwise', ...) {
  type <- c('pointwise')[min(pmatch(type,c('pointwise')))]

  if(type=='pointwise'){
    return( apply(object$data, MARGIN = 1, sd) )
  } else{
    stop('Only type="pointwise" is implemented', call. = FALSE)
  }

}


#' @rdname sdvar
#' @export
var <- function(object, ...) UseMethod("var")
#' @rdname sdvar
#' @export
var.default <- function(object, ...) stats::var(object, ...)


#' @rdname sdvar
#'
#' @export
#'
#' @examples
#' var(electricity,type='pointwise')
#' var(electricity,type='operator')
var.dfts <- function(object, type=c('operator','pointwise'), ...) {
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


#' @name dfts_group
#'
#' @return dfts object with data as cumsum
#' @export
#'
#' @examples
#' cumsum(electricity)
cumsum.dfts <- function(x){
  x <- dfts(x)

  dfts(t(apply(x$data,MARGIN = 1,cumsum)),
        labels = x$labels,
        fparam = x$fparam)
}
