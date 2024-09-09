#' Define funts Object
#'
#' @param x Data.frame to convert into data
#' @param labels Labels for the observations. Defaults to the column names of
#'  \code{x}.
#' @param intraobs Vector of numerics indicating the points for each
#'  observation. Defaults to even spacing on [0,1], but may be uneven or
#'  for any real numbers.
#'
#' @return funts object
#' @export
#'
#' @examples
#' funts(electricity)
#' funts(generate_brownian_motion(100, c(0,0.1,0.25,0.5,1)))
funts <- function(X, labels=colnames(as.data.frame(X)),
                  intraobs=seq(0,1,length.out=nrow(X))){
  if(is.funts(X)) return(X,labels=labels,intraobs=intraobs)
  # Note: as.matrix is very important! Otherwise extraction and computation is
  #   far more expensive!
  funts_obj <- list(
    'data' = as.matrix(X),
    'labels' = labels,
    'intraobs' = intraobs
  )
  colnames(funts_obj$data) <- NULL
  # structure(x, 'labels' = labels, class = "funts")
  class(funts_obj) <- 'funts'

  funts_obj
}


#' Check data to see if its funts
#'
#' @param x Data to examine class
#'
#' @return Boolean indicating if \code{x} is a funts object or not
#' @export
#'
#' @examples
#' is.funts(electricity)
#' is.funts(funts(electricity))
is.funts <- function(x){
  inherits(x, "funts") #& length(x$x) > 0
}


#' Convert to funts
#'
#' @param x Object to convert to funts
#'
#' @return funts object
#' @export
#'
#' @examples
#' as.funts(electricity)
as.funts <- function(x){
  switch(class(x)[[1]],
         'data.frame'={
           funts(x)
         },
         'matrix'={
           funts(x)
         },
         {
           stop('Only matrix and data.frame objects currently setup')
         }
  )
}


#' Validate funts
#'
#' @param x funts object to verify for any violations
#'
#' @return Silently returns the data if there are no violations
#' @export
#'
#' @examples
#' validate.funts(funts(electricity))
validate.funts <- function(x){
  if (!inherits(x, "funts")) stop("Objects must be of class funts")
  values <- unclass(x)$data
  labels <- x$labels

  if (length(labels) != ncol(values)) {
    stop("There must be the same number of `labels` as values in `x`")
  }

  invisible(x)
}


#' Math Group for funts
#'
#' @param x funts object
#' @param ... Additional parameters for functions in math group
#'
#' @return funts object with applied functions
#' @export
#'
#' @examples
#' sqrt( funts(electricity) )
Math.funts <- function (x, ...) {
  validate.funts(x)

  x$data <- callGeneric(x$data,...)

  x
}


#' Ops Group for funts
#'
#' TODO: Figure out for non-intersections.
#'
#' @param e1 funts object
#' @param e2 funts object
#'
#' @return funts object from the two elements
#' @export
#'
#' @examples
#' funts(electricity) + funts(electricity)
#' funts(electricity) * funts(electricity)
Ops.funts <- function (e1, e2) {
  # Combines the columns that are the same

  if(missing(e2)) {
    .Class <- "matrix"
    NextMethod()
  } else {
    c.e1 <- class(e1)
    c.e2 <- class(e2)

    if("funts" %in% c.e1 && "funts" %in% c.e2) {
      nce1 <- NCOL(e1$data)
      nce2 <- NCOL(e2$data)

      if(nce1 != nce2){#} && nce1 != 1 && nce2 != 1) {
        stop("Ops.funts: non conformable data.")
      }
      stopifnot( all.equal(class(e1$labels),
                           class(e2$labels)) )
      stopifnot( all.equal(class(e1$intraobs),
                           class(e2$intraobs)) )

      i.cols <- intersect( e1$labels, e2$labels )

      class(i.cols) <- class(e1$labels)

      ## if there is an intersection, the do the Op
      if(length(i.cols)) {
        io <- e1$intraobs
        e1 <- e1$data[,e1$labels==i.cols]
        e2 <- e2$data[,e2$labels==i.cols]
        # e1 <- e1[i.dates,]
        # e2 <- e2[i.dates,]

        # if(nce1==1 && nce2!=1) {
        #   e1 <- rep(e1,nce2)
        # } else if(nce1!=1 && nce2==1) {
        #   e2 <- rep(e2,nce1)
        # }

        .Class <- "matrix"
        ans <- NextMethod()

        #colnames(ans) <- i.cols
        ans <- funts(ans,labels = i.cols,intraobs = io)
      } else {
        ## no matching cols, return NULL
        ans <- NULL
      }
    } else {
      .Class <- "matrix"
      ans <- NextMethod()
      colnames(ans) <- i.cols
      ans <- funts(ans,labels = i.cols)

      # if("funts" %in% c.e1) {
      #   ans.dates <- attr(e1,"index")
      # } else {
      #   ans.dates <- attr(e2,"index")
      # }
      # attr(ans,"index") <- ans.dates
      # class(ans) <- c("funts","zoo")
    }
    ans
  }
}


#' Difference funts
#'
#' @param x funts object
#' @param lag An integer indicating which lag to use.
#' @param differences	An integer indicating the order of the difference.
#' @param ... Further arguments to be passed to or from methods.
#'
#' @return funts object with differenced values
#' @export
#'
#' @examples
#' diff(funts(electricity), lag=1)
#' diff(funts(electricity), differences=2)
diff.funts <- function(x, lag = 1L, differences = 1L, ...) {
  if(differences==0) return(x)

  # Data check
  if (length(lag) != 1L || length(differences) > 1L || lag <
      1L || differences < 1L)
    stop("'lag' and 'differences' must be integers >= 1")
  if (lag * differences >= ncol(x$data))
    return(NULL)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")

  dat1 <- data.frame(matrix(NA,
                            nrow = nrow(x$data),
                            ncol = ncol(x$data)-lag))
  #colnames(dat1) <- x$labels[-c(length(x$labels)-(lag-1):0)]
  for(i in 1:(ncol(x$data)-lag)){
    dat1[,i] <- x$data[,i+lag]-x$data[,i]
  }

  diff.funts(x=funts(dat1),lag=lag, differences=differences-1)
}


#' Lag funts
#'
#' @param x funts object
#' @param k integer indicating the number of lags for the data
#' @param ... Unused additional parameters
#'
#' @return A funts object
#' @export
#'
#' @examples
#' lag.funts(funts(electricity))
lag.funts <- function(x, k, ...) {
  stopifnot(k >= 0)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")
  #
  dat1 <- x$data[,-c(1:round(k))]
  funts(dat1,labels = x$labels[-c(1:round(k))])
}


#' Max funts
#'
#' Get the observation(s) with the max average value.
#'
#' TODO:: Add a pw option
#'
#' @param x funts object
#' @param ... Unused additional parameters
#'
#' @return A funts object
#' @export
#'
#' @examples
#' max(funts(electricity))
max.funts <- function(x, ...){
  idx <- which.max(colMeans(x$data))
  funts(x$data[,idx,drop=FALSE],
        labels = x$labels[idx])
}


#' Min funts
#'
#' Get the observation(s) with the min average value
#'
#' TODO:: Add a pw option
#'
#' @inheritParams max.funts
#'
#' @return A funts object
#' @export
#'
#' @examples
#' min(funts(electricity))
min.funts <- function(x, ...){
  idx <- which.min(colMeans(abs(x$data)))
  funts(x$data[,idx,drop=FALSE],
        labels = x$labels[idx])
}


#' Mean for funts
#'
#' @inheritParams max.funts
#'
#' @return Numeric for the mean at each resolution
#' @export
#'
#' @examples
#' mean(funts(electricity))
mean.funts <- function(x, ...) {
  rowMeans(x$data)
}


#' Median for funts
#'
#' @inheritParams max.funts
#'
#' @return Numeric for the median at each resolution
#' @export
#'
#' @examples
#' median(funts(electricity))
median.funts <- function(x, ...) {
  apply(x$data, MARGIN = 1, median)
}


#' Dimension of funts
#'
#' @inheritParams max.funts
#'
#' @return Numerics indicating the dimension of the funts object
#' @export
#'
#' @examples
#' dim(funts(electricity))
dim.funts <- function(x, ...) {
  dim(x$data)
}


#' Plot funts
#'
#' @param x A funts object
#' @param ... Additional parameters to pass into plotting
#'
#' @return Plot
#' @export
#'
#' @examples
#' plot(funts(electricity))
plot.funts <- function(x, ...){
  plot_fd(x$data, ...)
}


#' Quantile funts
#'
#' @param x A funts object
#' @param ... Additional parameters to pass into quantile function
#'
#' @return Matrix with columns for each requested quantile
#' @export
#'
#' @examples
#' quantile(funts(electricity))
#' quantile(funts(electricity),probs = 0.95)
quantile.funts <- function(x, probs = seq(0, 1, 0.25), ...){
  if(length(probs)>1){
    output <- t(apply(x$data, MARGIN = 1, quantile, probs=probs, ...))
  } else{
    output <- data.frame(apply(x$data, MARGIN = 1, quantile, probs=probs, ...))
    colnames(output) <- paste0(probs*100,'%')
  }

  output
}
