#' dfts Objects
#'
#' Define data as, check data if, and convert data to dfts object. The data
#'  may be from another popular implementation.
#'
#' @param X Data to convert into dfts object. Options include: data.frame,
#'  matrix, array, fda::fd, fda.usc::fdata, rainbow::fts (used in ftsa),
#'  rainbow::fds (used in ftsa), funData::funData, and fChange::dfts. For a
#'  matrix, each column is a unique observation, at the rows are the observed
#'  intra-observation (i.e. resolution) points.
#' @param name String for the name of the object. Defaults to the name of the
#'  input variable.
#' @param labels Labels for the observations. Defaults to the column names or
#'  names inside of the object X.
#' @param intratime Vector of numerics indicating the points for each
#'  observation. Defaults to even spacing on \[0,1\], or those included in the
#'  object. These may be unevenly spaced.
#' @param inc.warnings Boolean on if warnings should be given. Defaults to TRUE,
#'
#' @name dfts
#'
#' @return dfts / as.dfts: dfts object
#' @export
#'
#' @examples
#' bm <- dfts(generate_brownian_motion(100, c(0,0.1,0.25,0.5,1)))
#'
#' # Electricity
#' result <- dfts(electricity)
dfts <- function(X, name=NULL,
                  labels=NULL,
                  intratime=NULL,
                  inc.warnings=TRUE){
  ## Return if already right object
  if(is.dfts(X)) {
    if(inc.warnings & sum(is.na(X$data))>0)
      warning('NA values in data may affect some methods',call. = FALSE)

    return(.validate.dfts(X))
  }

  ## Fill-in Data
  if(is.null(name))
    name <- deparse(substitute(X))

  # Setup Data in Different Scenarios
  # Note: as.matrix is very important in data!
  #   Otherwise computation is far more expensive!
  if(requireNamespace("fda", quietly = TRUE) &&
     fda::is.fd(X)){
    if(is.null(labels))
      labels <- X$fdnames$reps
    if(is.null(intratime)){
      ran_vals <- X$basis$rangeval
      intratime <- seq(ran_vals[1],ran_vals[2],
                      length.out=length(X$fdnames$time))
    }

    dfts_obj <- list(
      'data' = as.matrix(fda::eval.fd(intratime, fdobj = X)),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else if(requireNamespace("fda.usc", quietly = TRUE) &&
            fda.usc::is.fdata(X)){
    if(is.null(labels))
      labels <- rownames(X$data)
    if(is.null(intratime))
      intratime <- X$argvals

    dfts_obj <- list(
      'data' = as.matrix(t(X$data)),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else if(requireNamespace('rainbow', quietly = TRUE) &&
            (inherits(X, "fts") || inherits(X, "fds"))){
    if(is.null(labels)){
      if(inherits(X, "fts"))
        labels <- X$time
      if(inherits(X, "fds"))
        labels <- colnames(X$y)
    }
    if(is.null(intratime))
      intratime <- X$x

    dfts_obj <- list(
      'data' = as.matrix(X$y),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else if(requireNamespace('funData', quietly = TRUE) &&
            inherits(X, "funData")){
    if(is.null(labels))
      labels <- rownames(X@X)
    if(is.null(intratime))
      intratime <- X@argvals[[1]]

    dfts_obj <- list(
      'data' = as.matrix(t(X@X)),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else if(inherits(X, "matrix") ||
            inherits(X, "data.frame") ||
            inherits(X, "array")){
    if(is.null(intratime))
      intratime <- seq(0,1,length.out=nrow(X))
    if(is.null(labels))
      labels <- colnames(as.data.frame(X))

    dfts_obj <- list(
      'data' = as.matrix(X),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else if(inherits(X, "numeric")){
    if(is.null(intratime))
      intratime <- seq(0,1,length.out=1)
    if(is.null(labels)){
      if(is.null(names(X))){
        labels <- 1:length(X)
      }else{
        labels <- names(X)
      }
    }
    dfts_obj <- list(
      'data' = matrix(X,nrow=1),
      'name' = name,
      'labels' = labels,
      'intratime' = intratime
    )
  } else{
    stop('Data type not supported.',call. = FALSE)
  }

  ## Clean-up
  #   Old: structure(x, 'labels' = labels, class = "dfts")
  colnames(dfts_obj$data) <- NULL
  rownames(dfts_obj$data) <- NULL
  class(dfts_obj) <- 'dfts'

  ## Verification
  # Note if NAs are present
  if(length(dfts_obj$intratime)!=nrow(dfts_obj$data))
    stop('Resolution does not match data (rows)',call. = FALSE)
  if(length(dfts_obj$labels)!=ncol(dfts_obj$data))
    warning('Labels do not match data (cols)',call. = FALSE)

  if(inc.warnings & sum(is.na(dfts_obj$data))>0)
    warning('NA values in data may affect some methods',call. = FALSE)

  .validate.dfts(dfts_obj)
}


#' @rdname dfts
#' @export
#'
#' @examples
#' result <- as.dfts(electricity)
as.dfts <- function(X){
  dfts(X)
}


#' @rdname dfts
#' @export
#'
#' @return is.dfts: Boolean indicating if \code{x} is a dfts object or not
#'
#' @examples
#' result <- is.dfts(electricity)
is.dfts <- function(X){
  inherits(X, "dfts")
}



#' Validate dfts
#'
#' Internal function to check lengths and data information is right
#'
#' @param X dfts object to verify for any violations
#'
#' @return Silently returns the data if there are no violations
#'
#' @keywords internal
#' @noRd
.validate.dfts <- function(X){
  # Class
  if (!inherits(X, "dfts")) stop("Objects must be of class dfts")

  # Sizes
  if(ncol(X) != length(X$labels))
    stop("Data must have the same number of columns as the length of labels")
  if(nrow(X) != length(X$intratime))
    stop("Data must have the same number of rows as the length of intratime")

  X
}



#' Group Generic Functions
#'
#' Group generic methods defined for things like Math, Ops, and so forth.
#'
#' @param x,e1,e2 dfts object
#' @param ... Further arguments passed to the methods
#'
#' @return dfts object
#' @export
#'
#' @name dfts_group
#'
#' @examples
#' result <- sqrt( electricity )
Math.dfts <- function (x, ...) {
  x$data <- methods::callGeneric(x$data,...)

  x
}


#' @rdname dfts_group
#' @export
#'
#' @examples
#' result <- electricity + electricity
#' result1 <- electricity * electricity
Ops.dfts <- function (e1, e2) {
  # Combines the columns that are the same
  # TODO: Figure out for non-intersections.

  if(missing(e2)) {
    .Class <- "matrix"
    NextMethod()
  } else {
    c.e1 <- class(e1)
    c.e2 <- class(e2)

    if("dfts" %in% c.e1 && "dfts" %in% c.e2) {
      nce1 <- NCOL(e1$data)
      nce2 <- NCOL(e2$data)

      if(nce1 != nce2){#} && nce1 != 1 && nce2 != 1) {
        stop("Ops.dfts: non conformable data.")
      }
      stopifnot( all.equal(class(e1$labels),
                           class(e2$labels)) )
      stopifnot( all.equal(class(e1$intratime),
                           class(e2$intratime)) )

      i.cols <- intersect( e1$labels, e2$labels )

      class(i.cols) <- class(e1$labels)

      ## if there is an intersection, the do the Op
      if(length(i.cols)) {
        io <- e1$intratime
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
        ans <- dfts(ans,labels = i.cols,intratime = io)
      } else {
        ## no matching cols, return NULL
        ans <- NULL
      }
    } else {
      .Class <- "matrix"
      ans <- NextMethod()
      colnames(ans) <- i.cols
      ans <- dfts(ans,labels = i.cols)

      # if("dfts" %in% c.e1) {
      #   ans.dates <- attr(e1,"index")
      # } else {
      #   ans.dates <- attr(e2,"index")
      # }
      # attr(ans,"index") <- ans.dates
      # class(ans) <- c("dfts","zoo")
    }
    ans
  }
}


#' Difference dfts
#'
#' Difference the functional data at some lag / to some degree.
#'
#' @param x A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param lag An integer indicating which lag to use.
#' @param differences	An integer indicating the order of the difference.
#' @param ... Further arguments to be passed to or from methods.
#'
#' @return dfts object with differenced values
#' @export
#'
#' @examples
#' result <- diff(electricity, lag=1)
#' result1 <- diff(electricity, differences=2)
diff.dfts <- function(x, lag = 1L, differences = 1L, ...) {
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
  for(i in 1:(ncol(x$data)-lag)){
    dat1[,i] <- x$data[,i+lag]-x$data[,i]
  }

  diff.dfts(x=dfts(dat1),lag=lag, differences=differences-1)
}


#' Max / Min for dfts Objects
#'
#' Get the observation(s) with the min / max values. Selected as the observation
#'  with the largest / smallest mean or the pointwise values.
#'
#' @param x A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param type String indicating if finding for observation ('Obs', default),
#'  or for pointwise values ('intratime').
#' @param na.rm Boolean if NA values should be removed. Defaults to TRUE.
#' @param ... Additional parameters to stats function
#'
#' @return A dfts object
#' @export
#'
#' @name minmax
#'
#' @examples
#' results <- max(electricity)
max.dfts <- function(x, type=c('Obs','intratime'), na.rm=TRUE, ...){
  poss_types <- c('Obs','intratime')
  type <- .verify_input(type, poss_types)
  x <- dfts(x)

  if(type == 'Obs'){
    idx <- which.max(colMeans(x$data,na.rm = na.rm))
    return_dat <- dfts(x$data[,idx,drop=FALSE],
                        name = x$name,
                        intratime = x$intratime,
                        labels = x$labels[idx])
  }else if(type=='intratime'){
    dat <- apply(x$data,MARGIN = 1, max, na.rm=na.rm, ...)
    return_dat <- dfts(matrix(dat),
                        name = x$name,
                        intratime = x$intratime)

  }

  return_dat
}


#' @rdname minmax
#' @export
#'
#' @examples
#' results <- min(electricity)
min.dfts <- function(x, type=c('Obs','intratime'), na.rm=TRUE, ...){
  poss_types <- c('Obs','intratime')
  type <- .verify_input(type, poss_types)

  if(type == 'Obs'){
    idx <- which.min(colMeans(x$data,na.rm=na.rm))
    return_dat <- dfts(x$data[,idx,drop=FALSE],
                        name = x$name,
                        intratime = x$intratime,
                        labels = x$labels[idx])
  }else if(type=='intratime'){
    dat <- apply(x$data,MARGIN = 1, min, na.rm=na.rm, ...)
    return_dat <- dfts(matrix(dat),
                        name = x$name,
                        intratime = x$intratime)

  }

  return_dat
}


#' Average Functions for dfts Objects
#'
#' Compute the average values for dfts objects
#'
#' @inheritParams minmax
#'
#' @return Numeric vector
#' @export
#'
#' @name average
#'
#' @examples
#' results <- mean(electricity)
mean.dfts <- function(x, na.rm=TRUE, ...) {
  rowMeans(x$data,na.rm = na.rm, ...)
}


#' @rdname average
#' @export
#'
#' @importFrom stats median
#'
#' @examples
#' results <- median(electricity)
median.dfts <- function(x, na.rm=TRUE, ...) {
  apply(x$data,MARGIN = 1, stats::median, na.rm=na.rm, ...)
}


#' Quantile dfts
#'
#' Obtain the pointwise quantile information of dfts objects
#'
#' @param x A dfts object
#' @param probs Numerics in \[0,1\] indicating the probabilities of interest
#' @param ... Additional parameters to pass into quantile function
#'
#' @return Matrix with columns for each requested quantile
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' result <- quantile(electricity)
#' result1 <- quantile(electricity,probs = 0.95)
quantile.dfts <- function(x, probs = seq(0, 1, 0.25), ...){
  if(length(probs)>1){
    output <- t(apply(x$data, MARGIN = 1, stats::quantile, probs=probs, ...))
  } else{
    output <- data.frame(apply(x$data, MARGIN = 1, stats::quantile, probs=probs, ...))
    colnames(output) <- paste0(probs*100,'%')
  }

  output
}


#' Lag dfts
#'
#' @param x A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param k integer indicating the number of lags for the data
#' @param ... Unused additional parameters
#'
#' @return A dfts object
#' @export
#'
#' @importFrom stats lag
#'
#' @examples
#' result <- lag(electricity)
lag.dfts <- function(x, k=1, ...) {
  stopifnot(k >= 0)
  #ans <- .Call("lag", x, as.integer(k),PACKAGE="fts")
  x <- dfts(x)

  dat1 <- x$data[,-c(1:round(k))]
  dfts(dat1, labels = x$labels[-c(1:round(k))])
}


#' Dimension of dfts Object
#'
#' Retrieve the dimension of a dfts object.
#'
#' @inheritParams minmax
#'
#' @return Numerics indicating the dimension of the dfts object
#' @export
#'
#' @examples
#' results <- dim(electricity)
dim.dfts <- function(x, ...) {
  dim(x$data)
}

