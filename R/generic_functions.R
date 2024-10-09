#' Generic ACF/PACF
#'
#' @param x Object for computation of (partial) autocorrelation function
#'  (ACF/PACF)
#' @param ... Additional parameters to function based on data
#'
#' @name acf
#'
#' @return ACF or PACF values and plots
#' @export
#'
#' @examples
#' acf(1:10)
NULL

#' Generic Centering of Data
#'
#' @param object Object for computation of centering
#' @param ... Additional parameters to function based on data
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

#' @rdname acf
#'
#' @export
acf <- function(x, ...) UseMethod("acf")
#' @rdname acf
#'
#' @export
acf.default <- function(x, ...) stats::acf(x)


#' @rdname acf
#'
#' @export
pacf <- function(x, ...) UseMethod("pacf")
#' @rdname acf
#'
#' @export
pacf.default <- function(x, ...) stats::pacf(x)


#' ACF for Functional Data
#'
#' @inheritParams acf
#'
#' @return Functional ACF values
#' @export
#'
#' @examples
#' acf(funts(electricity))
acf.funts <- function(x, ...){
  invisible(.compute_FACF(x, ...))
}


#' PACF for Functional Data
#'
#' @inheritParams pacf
#'
#' @return Functional PACF values
#' @export
#'
#' @examples
#' pacf(funts(electricity))
pacf.funts <- function(x, ...){
  invisible(.compute_FPACF(x, ...))
    # n_harm = NULL, lag.max = NULL, ci=0.95, figure = TRUE, ...)
}


#' @rdname center
#' @export
center <- function(x, ...) UseMethod("center")
#' @rdname center
#' @export
center.default <- function(x, ...) { x - mean(x) }
#' @rdname center
#' @export
center.data.frame <- function(x, ...) { x <- as.matrix(x); NextMethod("center") }
#' @rdname center
#' @export
center.matrix <- function(x, ...) { x - rowMeans(x) }

#' Center funts data
#'
#' @rdname center
#'
#' @export
center.funts <- function(x, CPs=NULL, type='mean', ...) {
  type <- .verify_input(type, c('mean','median'))
  if(type=='mean'){
    row_method <- rowMeans
  }else if(type=='median'){
    row_method <- function(y){apply(y, MARGIN = 1, median)}
  }


  if(is.null(CPs)){
    x$data <- x$data - row_method(x$data)
  }else{
    CPs_use <- unique(c(0, CPs, ncol(x$data)))
    CPs_use <- CPs_use[order(CPs_use)]
    for(i in 1:(length(CPs)-1)){
      x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] <-
        x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] -
        row_method(x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE])
    }
  }
  x
}



#' @rdname pca
#' @export
pca <- function(x, ...) UseMethod("pca")
#' @rdname pca
#' @export
pca.default <- function(x, ...) { prcomp(x, ...) }


#' Principal component analysis
#'
#' @inheritParams pca
#' @param TVE Numeric in [0.1] for the total variance explained. Can use this to
#'  select only the required components
#' @param ... Additional parameters for \code{prcomp}.
#'
#' @return Principal component data
#' @export
#'
#' @examples
#' pca(funts(electricity))
pca.funts <- function(x, TVE = 1, ...){
  if(TVE > 1 || TVE < 0) stop('TVE must be in [0,1].',call. = FALSE)

  pc <- stats::prcomp(x=t(x$data), ...)
  if(TVE==1){
    min_pc <- length(pc$sdev)
  } else{
    min_pc <- min( which(cumsum(pc$sdev^2)/sum(pc$sdev^2) > TVE) )
  }

  ## TODO:: Check out
  new_rot <- pc$rotation[,1:min_pc,drop=FALSE]*sqrt(nrow(x$data))

  # scores <- matrix(NA, nrow=ncol(x$data),ncol = min_pc)
  # for(i in 1:ncol(x$data)){
  #   for(l in 1:min_pc){
  #     scores[i,l] <- dot_integrate(x$data[,i] %*% t(new_rot[,l]))
  #   }
  # }
  scores <- (t(x$data-pc$center) %*% new_rot)[,1:min_pc,drop=FALSE]/nrow(x$data)

  # res = fda::pca.fd(fda::Data2fd(X$data))
  # pc$sdev^2
  # res$values[1:2]
  # View(pc$x); View(res$scores)
  # View(pc$rotation); View(fda::eval.fd(seq(0,1,length.out=30),res$harmonics))

  list(sdev = pc$sdev[1:min_pc]/sqrt(nrow(x$data)),
       rotation = new_rot,
       center = pc$center,
       scale = pc$scale,
       x = scores,
       orig_pc = pc)
}


#' @rdname sd
#' @export
sd <- function(x, ...) UseMethod("sd")
#' @rdname sd
#' @export
sd.default <- function(x, ...) stats::sd(x, ...)

#' SD for funts
#'
#' @inheritParams sd
#' @param type Character to specify if an operator ('op') or pointwise ('pw')
#'  calculation is desired
#' @param ... Unused parameters
#'
#' @return
#' @export
#'
#' @examples
#' sd(funts(electricity),type='pointwise')
sd.funts <- function(x, type='pointwise', ...) {
  type <- c('pointwise')[min(pmatch(type,c('pointwise')))]

  if(type=='pointwise'){
    return( apply(x$data, MARGIN = 1, sd) )
  } else{
    stop('Only type="pointwise" is implemented', call. = FALSE)
  }

}


#' @rdname sd
#' @export
var <- function(x, ...) UseMethod("var")
#' @rdname sd
#' @export
var.default <- function(x, ...) stats::var(x, ...)

#' Variance for funts
#'
#' @inheritParams sd.funts
#'
#' @return
#' @export
#'
#' @examples
#' var(funts(electricity),type='pointwise')
var.funts <- function(x, type=c('operator','pointwise'), ...) {
  type <- c('operator','pointwise')[min(pmatch(type,c('operator','pointwise')))]

  if(type=='operator'){
    autocov_approx_h(X$data,0)
  } else if(type=='pointwise'){
    Apply(x$data, MARGIN = 1, var)
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
cumsum.funts <- function(x){
    funts(t(apply(x$data,MARGIN = 1,cumsum)),
          labels = x$labels,
          intraobs = x$intraobs)
}
