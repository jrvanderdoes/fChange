#' Define funts Object
#'
#' @param x Data.frame to convert into data
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
#' funts(electricity)
funts <- function(X, labels=colnames(as.data.frame(X)),
                  intraobs=seq(0,1,length.out=nrow(X))){
  funts_obj <- list(
    'data' = as.data.frame(X),
    'labels' = labels,
    'intraobs' = intraobs
  )
  #structure(x, 'labels' = labels, class = "funts")
  class(funts_obj) <- 'funts'

  funts_obj
}


#' Check is funts
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
is.funts <- function(x){
  inherits(x, "funts") #& length(x$x) > 0
}


#' Check as funts
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
#' @param x funts to verify class does not have violations
#'
#' @return
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

  # return(x)
}


#' Math Group for funts
#'
#' @param x
#' @param ...
#'
#' @return
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
#' @param e1
#' @param e2
#'
#' @return
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

      i.cols <- intersect( e1$labels, e2$labels )

      class(i.cols) <- class(e1$labels)

      ## if there is an intersection, the do the Op
      if(length(i.cols)) {
        e1 <- e1$data[,i.cols]
        e2 <- e2$data[,i.cols]
        # e1 <- e1[i.dates,]
        # e2 <- e2[i.dates,]

        # if(nce1==1 && nce2!=1) {
        #   e1 <- rep(e1,nce2)
        # } else if(nce1!=1 && nce2==1) {
        #   e2 <- rep(e2,nce1)
        # }

        .Class <- "data.frame"
        ans <- NextMethod()
        colnames(ans) <- i.cols
        funts(ans,labels = i.cols)
      } else {
        ## no matching cols, return NULL
        ans <- NULL
      }
    } else {
      .Class <- "data.frame"
      ans <- NextMethod()
      colnames(ans) <- i.cols
      funts(ans,labels = i.cols)

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
#' @param x
#' @param k
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' diff(funts(electricity), k=1)
#' diff(funts(electricity), k=2)
diff.funts <- function(x, k, ...) {
  stopifnot(k > -1)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")

  dat1 <- data.frame(matrix(NA,
                            nrow = nrow(x$data),
                            ncol = ncol(x$data)-k))
  colnames(dat1) <- x$labels[-c(length(x$labels)-(k-1):0)]
  for(i in 1:(ncol(x$data)-k)){
    dat1[,i] <- x$data[,i+k]-x$data[,i]
  }

  funts(dat1)
}


#' Lag funts
#'
#' @param x
#' @param k
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
lag.funts <- function(x, k, ...) {
  stopifnot(k > -1)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")
  #
  dat1 <- x$data[,-c(1:k)]
  funts(dat1,labels = x$labels[-c(1:k)])
}


#' Max funts
#'
#' Get observations with the max absolute value
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' max(funts(electricity))
max.funts <- function(x, ...){
  idx <- which.max(colMeans(abs(x$data)))
  funts(x$data[,idx,drop=FALSE],
        labels = x$labels[idx])
}


#' Min funts
#'
#' Get observations with the min absolute value
#'
#' @param x
#' @param ...
#'
#' @return
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
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' mean(funts(electricity))
mean.funts <- function(x, ...) {
  rowMeans(x$data)
}


#' Median for funts
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' median(funts(electricity))
median.funts <- function(x, ...) {
  apply(x$data, MARGIN = 1, median)
}


#' SD for funts
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' sd(funts(electricity))
sd.funts <- function(x, ...) {
  apply(x$data, MARGIN = 1, sd)
}


#' Variance for funts
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' var(funts(electricity))
sd.funts <- function(x, ...) {
  apply(x$data, MARGIN = 1, var)
}


#' Dimension of funts
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' dim(funts(electricity))
dim.funts <- function(x, ...) {
  dim(x$data)
}


#' Plot funts
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' plot(funts(electricity))
plot.funts <- function(x, ...){
  plot_fd(x$data, ...)
}


#' Quantile funts
#'
#' @param x
#' @param ...
#'
#' @return
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
