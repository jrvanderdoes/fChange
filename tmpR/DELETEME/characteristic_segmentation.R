#' Complete Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. Change points are recursively found until no
#'     more change points are detected.
#'
#' @param X funts object or numeric data.frame with rows for evaluated values and columns
#'    indicating functional observations
#' @param statistic Type of test statistic, \code{Tn} (integrated) or \code{Mn} (maximized)
#' @param method description
#' @param M description
#' @param J description
#' @param space description
#' @param h_function description
#' @param iters description
#' @param replace description
#' @param alpha Numeric value in \eqn{[0, 1]} indicating the significance for
#'     cutoff_function.
#' @param final_verify (Optional) Boolean value XXXXXXXXXXXXX
#'
#' Indicates if a final pass looking at sequences with only one change point
#'     should be conducted to verify results. Note, this may modify existing
#'     locations of change points, potentially to less accurate locations.
#' @param silent (Optional) Boolean value
#'
#' Indicates if useful output should be silenced. Default FALSE shows output.
#'
#' @return A list of numeric values indicating change points  (if exists),
#'     NA otherwise
#' @export
#'
#' @examples
#' # binary_segmentation(funts(electricity[,1:50]),"Tn","Boot")
binary_segmentation <- function(X, statistic=c('Tn','Mn'),
                                method=c('Sim','Approx','Boot'),
                                M=20, J=50, space='BM',
                                h_function = function(X){ncol(X)^(1/3)},
                                iters = 1000,
                                replace = FALSE, alpha = 0.05,
                                final_verify = TRUE,
                                silent = FALSE) {
  X <- funts(X)

  # Test Statistics
  if(length(statistic)!=1 || !(statistic %in% c('Tn','Mn'))){
    stop('Choose "Tn" or "Mn" as the test statistic')
  }else if(length(method)!=1 || !(method %in% c('Sim','Approx','Boot'))){
    stop('Choose "Sim", "Approx", or "Boot" as the method')
  }

  if(statistic=='Tn'){
    if(method=='Sim'){

      .fn = function(X,
                    M, J, space,
                    blocksize, iters,
                    alpha, h_function,
                    replace,
                    K = bartlett_kernel,
                    ...){
        X <- funts(X)

        h <- h_function(X$data)

        W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

        tmp <-.ce_detect_Tn(X=X$data, W=W, M=M, J=J, nSims=iters,
                            h=h, K=bartlett_kernel, silent=TRUE)

        # if(return_pval){
        #   cp <- ifelse(tmp$pval<=alpha, compute_Mn(X$data,W=W,J=J)$location, NA)
        #   return(c(cp,tmp$pval))
        # }

        ifelse(tmp$pval<=alpha, compute_Mn(X$data,W=W,J=J)$location, NA)
      }

    }else if(method=='Approx'){

      .fn = function(X,
                    M, J, space,
                    blocksize, iters,
                    alpha, h_function,
                    replace,
                    K = bartlett_kernel,
                    ...){
        X <- funts(X)

        h <- h_function(X$data)

        W = computeSpaceMeasuringVectors(M=M,space=space,X=X)

        stat <- compute_Tn(X=X$data, W=W, J=J)
        cutoff <- .compute_Welch(X$data,alpha = alpha,W=W,M=M,h=h,K=K)

        ifelse(stat>=cutoff,compute_Mn(X$data,W=W,J=J)$location,NA)
      }

    }else{
      .fn = function(X,
                    M, J, space,
                    blocksize, iters,
                    alpha, h_function,
                    replace,
                    K = bartlett_kernel,
                    ...){
        X <- funts(X)

        h <- h_function(X$data)

        W = computeSpaceMeasuringVectors(M=M,space=space,X=X)

        ifelse(.ce_bootstrap(X=X$data,statistic='Tn',W=W,J=J,space=space,
                             blocksize=h,iters=iters,
                                        replace=replace,alpha=alpha,
                                        silent=TRUE)$pval <= alpha,
               compute_Mn(X$data,W=W,J=J)$location,NA)
      }

    }
  }else{
    if(method=='Sim'){

      .fn = function(X,
                    M, J, space,
                    blocksize, iters,
                    alpha, h_function,
                    replace,
                    K = bartlett_kernel,
                    ...){
        X <- funts(X)

        h <- h_function(X$data)

        W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

        tmp <-.ce_detect_Mn(X=X$data, W=W, M=M, J=J, nSims=iters,
                                           h=h, K=bartlett_kernel, silent=TRUE)

        ifelse(tmp$pval<=alpha,tmp$location,NA)
      }


    }else if(method=='Approx'){
      stop('Error no Approx method for Mn')
    }else{
      .fn = function(X,
                    M, J, space,
                    blocksize, iters,
                    alpha, h_function,
                    replace,
                    K = bartlett_kernel,
                    ...){
        X <- funts(X)

        h <- h_function(X$data)

        W = computeSpaceMeasuringVectors(M=M,space=space,X=X)

        ifelse(.ce_bootstrap(X=X$data,statistic='Mn',W=W,J=J,space=space,
                             blocksize=h,iters=iters,
                                             replace=replace,alpha=alpha,
                                             silent=TRUE)$pval <= alpha,
               compute_Mn(X$data,W=W,J=J)$location,NA)
      }

    }
  }

  # Get change points
  # TODO:: THis needs fixing
  blocksize <- K <- bartlett_kernel
  CPsVals <- .ce_recursive_segmentation(
    X=X$data, fn=.fn,
    M=M, J=J, space=space,
    blocksize=blocksize, K=K,
    h_function = h_function, iters = iters,
    replace = replace, alpha = alpha,
    silent = silent,
    addAmt = 0
  )

  # Verify as desired
  if (final_verify) {
    CPsVals <- .ce_verify_changes(
      CPsVals = CPsVals, X = X,
      fn=.fn, M=M, J=J, space=space,
      blocksize=blocksize, K=K,
      h_function=h_function, iters=iters,
      replace=replace,
      alpha=alpha,
      silent = silent
    )
  }

  return(CPsVals)
}


#' Characteristic Change - Recursive Segmentation
#'
#' @inheritParams binary_segmentation
#'
#' @param fn Method for test statistic
#' @param addAmt Number to add when recursively splitting to get correct
#'  change location
#' @param blocksize unused
#' @param K unused
#'
#' @return Numeric vector of changes
.ce_recursive_segmentation <- function(X, fn,
                        M=20, J=50, space='BM',
                        h_function = function(X){ncol(X)^(1/3)},
                        iters = 1000,
                        replace = FALSE, alpha = 0.05,
                        silent = FALSE,
                        addAmt = 0,
                        blocksize=NA, K=NA) {
  if(ncol(X)<=1) return(NA)

  h <- h_function(X)
  # Look for a single change
  potential_cp <- fn(X, M=M, J=J,
                     space=space,
                     h=h, iters=iters,
                     replace=replace, alpha=alpha
                    )

  # No Change Point Detected
  if (is.na(potential_cp)) {
    return()
  }

  # Display progress
  if (!silent) {
    cat(paste0(
      "ChangePoint Detected (", 1 + addAmt, "-", addAmt + ncol(X), " at ",
      addAmt + potential_cp, "): Segment Data and Re-Search\n"
    ))
  }

  # Search Recursively
  return(c(
    .ce_recursive_segmentation(
      X = as.data.frame(X[, 1:potential_cp]),
      fn=fn,
      #trim_function,
      M=M, J=J, space=space,
      h_function = h_function, iters = iters,
      replace = replace, alpha = alpha,
      silent = silent,
      addAmt = addAmt),
    potential_cp + addAmt,
    .ce_recursive_segmentation(
      X = as.data.frame(X[, (potential_cp + 1):ncol(X)]),
      fn=fn,
      #trim_function,
      M=M, J=J, space=space,
      h_function = h_function, iters = iters,
      replace = replace, alpha = alpha,
      silent = silent,
      addAmt = addAmt + potential_cp)
  ))
}


#' Characteristic Change - Verify Segmentation
#'
#' @inheritParams binary_segmentation
#' @param CPsVals Numeric vector of potential change points
#' @param fn Method for test statistic
#' @param ... Additional arguments passed into test_statistic_function
#'
#' @return Vector of change points or NA if none
#'
#' @keywords internal
#' @noRd
.ce_verify_changes <- function(CPsVals, X,
                               fn, M=M, J=J, space=space,
                               h_function=function(X){ncol(X)^(1/3)}, iters=iters,
                               replace=replace,
                               alpha=alpha,
                               silent = FALSE,
                               ...) {
  X <- funts(X)

  if (!silent) cat("-- Verify Step --\n")

  if (length(CPsVals) >= 1) { # If there was a CP detected
    tmp_cps <- c(0, CPsVals, ncol(X$data))
    tmp_cps <- tmp_cps[order(tmp_cps)]
    CPsVals <- pvals <- c()
    for (i in 2:(length(tmp_cps) - 1)) {
      # Get CP
      h <- h_function(X$data[,(tmp_cps[i-1]+1):tmp_cps[i+1]])
      potential_cp <- fn(X$data[,(tmp_cps[i-1]+1):tmp_cps[i+1]], M=M, J=J,
                         space=space,
                         h=h, iters=iters,
                         replace=replace, alpha=alpha,
                         return_pval=TRUE
      )

      if (!is.na(potential_cp[1])) {
        CPsVals <- c(CPsVals, potential_cp[1] + tmp_cps[i - 1])
        pvals <- c(pvals,potential_cp[2])
      }
    }
    if(length(CPsVals)>=1){
      return(list(
        CPsVals[order(CPsVals)],
        pvals[order(CPsVals)]))
    }else{
      return(NA)
    }
  } else {
    # Get CP
    h <- h_function(X$data)
    CPsVals <- fn(X$data, M=M, J=J,
                  space=space,
                  h=h, iters=iters,
                  replace=replace, alpha=alpha
    )
  }

  # Order and return
  CPsVals <- unlist(CPsVals)
  if (sum(is.na(CPsVals)) == length(CPsVals)) {
    return(NA)
  }
  if( length(CPsVals)<=1) return(CPsVals)

  CPsVals[order(CPsVals)]
}


#' Characteristic Detection - Tn Method
#'
#' @inheritParams binary_segmentation
#' @param W Space Measuring Vectors
#' @param nSims iters from before
#' @param h Window for kernel function
#' @param K Kernel function
#'
#' @return List with p-value for change detection
.ce_detect_Tn <- function(X, W,
                          M = 20, J=50,
                          nSims = 1000,
                          h = 3,
                          K = bartlett_kernel,
                          silent = FALSE) {
  X <- funts(X)

  # Determine Number of Iterations
  val_Tn <- compute_Tn(X$data, W, J)

  MJ <- M * J

  sqrtMat <- .compute_sqrtMat(X$data,W,J,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    # todo / ncol(mvnorms1)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    #apply(abs(gamVals)^2, MARGIN = 1, dot_integrate)
    apply(results, MARGIN = 1, dot_integrate)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Tn)
  )
}


#' Characteristic Detection - Mn Method
#'
#' @inheritParams .ce_detect_Tn
#'
#' @return List with p-value for change detection
.ce_detect_Mn <- function(X, W,
                          M = 20, J=50,
                          nSims = 1000,
                          h = 3,
                          K = bartlett_kernel,
                          silent = FALSE) {

  # Determine Number of Iterations
  val_Mn <- compute_Mn(X, W, J)

  MJ <- M * J

  sqrtMat <- .compute_sqrtMat(X,W,J,h,K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms),sqrtMat)
    # 0 0 0 0 (M times) ... .... 1 1 1 1 (M times)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    apply(results, MARGIN = 1, max)
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcess <- as.vector(gamProcess)

  list(
    'pval'=1 - stats::ecdf(gamProcess)(val_Mn$value),
    'location'=val_Mn$location
  )
}


#' Characteristic Detection - Mn Method
#'
#' @inheritParams .ce_recursive_segmentation
#' @inheritParams .ce_detect_Tn
#' @param statistic Test statistic of interest
#'
#' @return List with (1) value, (2) cutoff, (3) p-value, and (4) bootstrapped samples
.ce_bootstrap <- function(X, statistic=c('Tn','Mn'),
                          W=computeSpaceMeasuringVectors(M = 20, space=space, X=X),
                          J=50, space='BM',
                          blocksize=1, iters = 1000,type='separate',
                          replace = FALSE, alpha = 0.05,
                          silent = FALSE) {
  # Test Statistics
  if(length(statistic)!=1){
    stop('Choose "Tn" or "Mn" as the test statistic')
  }else if(statistic=='Tn'){
    fn <- compute_Tn
  }else if(statistic=='Mn'){
    fn <- compute_Mn
  }else{
    stop('Choose "Tn" or "Mn" as the test statistic')
  }

  ## Get Function Value and estimate time
  st <- Sys.time()
  full_val <- fn(X, W=W, J=J)[[1]]
  en <- Sys.time()

  if (!silent) {
    cat(paste0(
      "Estimated time: ",
      round(difftime(en, st, "units" = "mins")[[1]] * iters, 2),
      " mins\n"
    ))
  }

  # ## Create Permuted Samples
  # n <- ncol(X)
  # idxGroups <- .getChunks(1:n, n / blocksize)
  # idxs <- sapply(1:iters, function(i, m, indxs, replace) {
  #   samps <- sample(1:m, replace = replace)
  #   unlist(indxs[samps], use.names = FALSE)
  # }, m = length(idxGroups), indxs = idxGroups, replace = replace)
  # # If it is a matrix/dataframe this already even rows
  # #     (will be with permutation test)
  # if (!(is.matrix(idxs) | is.data.frame(idxs))) {
  #   idxs <- .convertSamplesToDF(idxs)
  # }
  #
  # ## Sample via bootstrap
  # bssamples <- sapply(as.data.frame(idxs),
  #                     function(loop_iter, fn, X1, W, J) {
  #                       fn(X1[, stats::na.omit(loop_iter)], W=W, J=J)[[1]]
  #                     },
  #                     X1 = X, fn = fn, W=W, J=J
  # )

  ## Create Permuted Samples
  bssamples <- .bootstrap(X,blocksize,iters,type,replace, fn=fn, W=W, J=J)

  list(
    "value" = full_val,
    "cutoff" = stats::quantile(bssamples, probs = c(1 - alpha))[[1]],
    "pval" = 1 - stats::ecdf(bssamples)(full_val),
    "BSSamples" = as.numeric(bssamples)
  )
}


