#' Generic ACF/PACF
#'
#' @param object Object for computation of (partial) autocorrelation function
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
#' @export
acf <- function(x, ...) UseMethod("acf")
#' @rdname acf
#' @inheritParams stats::acf
#' @export
acf.default <- function(x, ...) stats::acf(x)


#' @rdname acf
#' @export
pacf <- function(x, ...) UseMethod("pacf")
#' @rdname acf
#' @inheritParams stats::pacf
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
center.funts <- function(x, CPs=NULL, ...) {
  if(is.null(CPs)){
    x$data <- x$data - rowMeans(x$data)
  }else{
    CPs_use <- unique(c(0, CPs, ncol(x$data)))
    CPs_use <- CPs_use[order(CPs_use)]
    for(i in 1:(length(CPs)-1)){
      x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] <-
        x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE] -
        rowMeans(x$data[,(CPs_use[i]+1):CPs_use[i+1],drop=FALSE])
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

  if(TVE==1){
    return(stats::prcomp(x$data, ...))
  } else{
    pc <- stats::prcomp(x$data, ...)

    min_pc <- min(which(cumsum(pc$sdev)/sum(pc$sdev)>TVE))

    ## TODO:: Do I need sqrt(n)
    pc_new <- list(sdev = pc$sdev[1:min_pc]*sqrt(ncol(x$data)),
                   rotation = pc$rotation[,1:min_pc],#*sqrt(n),
                   center = pc$center,
                   scale = pc$scale,
                   x = pc$x[,1:min_pc],
                   fullpc = pc)
  }

  pc_new
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
#' sd(funts(electricity),type='pw')
sd.funts <- function(x, type=c('op','pw'), ...) {
  type <- c('op','pw')[min(pmatch(type,c('op','pw')))]

  if(type=='op'){
    stop('Not yet implemented, put type="pw" for implemented version',
         call. = FALSE)
  } else if(type=='pw'){
    return( apply(x$data, MARGIN = 1, sd) )
  }

  stop('Type must be "op" or "pw"', call. = FALSE)
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
#' var(funts(electricity),type='pw')
var.funts <- function(x, type=c('op','pw'), ...) {
  type <- c('op','pw')[min(pmatch(type,c('op','pw')))]

  if(type=='op'){
    stop('Not yet implemented, put type="pw" for implemented version',
         call. = FALSE)
  } else if(type=='pw'){
    return( apply(x$data, MARGIN = 1, var) )
  }

  stop('Type must be "op" or "pw"', call. = FALSE)
}
