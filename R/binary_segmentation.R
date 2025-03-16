#' Complete Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. Change points are recursively found until no
#'     more change points are detected.
#'
#' @inheritParams change
#' @param ... Additional arguments passed into method
#'
#' @return A data.frame of numeric values indicating change points and pvalues
#'
#' @noRd
#' @keywords internal
.binary_segmentation <- function(X, method,
                                 trim_function = function(X) {
                                   max(10, floor(log(ncol(X))),
                                       na.rm = TRUE
                                   )
                                 },
                                 alpha = 0.05, silent = FALSE, ...) {
  ## Setup
  X <- dfts(X)

  ## Get change points
  changes_info <- .detect_changes(
    X = X, method = method,
    trim_function = trim_function,
    alpha = alpha,
    addAmt = 0,
    silent = silent,
    ...
  )

  ## Verify
  changes1 <- .binary_verification(
    changes_info = changes_info, X = X, method = method,
    trim_function = trim_function, alpha = alpha,
    silent = silent, ...
  )

  ## Return Results
  changes_info
}



#' Detect Change Points
#'
#' This (internal) function is multiple .single_binary_segmentation for
#'     complete_binary_segmentation. It recursively calls itself
#'
#' @inheritParams change
#' @param ... Additional arguments passed into method
#'
#' @return A data.frame of numeric values indicating change points and pvalues
#'
#' @noRd
#' @keywords internal
.detect_changes <- function(X, method,
                            trim_function,
                            alpha = 0.05,
                            addAmt = 0,
                            silent = FALSE,
                            ...) {
  # Look for a single change
  potential <- .single_segment(X = X, method = method, trim_function = trim_function, ... )

  # No Change Point Detected
  return_now <- tryCatch({
    rval <- FALSE
    if (potential$pvalue > alpha) rval <- TRUE

    rval
  }, error = function(e){
    if (is.na(potential) || is.null(potential)) {
      TRUE
    }
  })
  if(return_now) return()

  # Display progress
  if (!silent) {
    cat(paste0(
      "ChangePoint Detected (", 1 + addAmt, "-", addAmt + ncol(X), " at ",
      addAmt + potential$location, "): Segment Data and Re-Search\n"
    ))
  }

  # Search Recursively
  return(
    rbind(
    .detect_changes(
      X = dfts(X$data[, 1:potential$location,drop=FALSE],
                labels = X$labels[1:potential$location],
                fparam = X$fparam,inc.warnings = FALSE),
      method = method,
      trim_function = trim_function,
      alpha = alpha,
      addAmt = addAmt,
      silent = silent,
      ...
    ),
    data.frame('location'=potential$location + addAmt, 'pvalue'=potential$pvalue),
    .detect_changes(
      X = dfts(X$data[, (potential$location + 1):ncol(X),drop=FALSE],
               labels = X$labels[(potential$location + 1):ncol(X)],
               fparam = X$fparam,inc.warnings = FALSE),
      method = method,
      trim_function = trim_function,
      alpha = alpha,
      addAmt = addAmt + potential$location,
      silent = silent,
      ...
    )
  ))
}


#' Single Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. At most one change point is detected.
#'
#'
#' @inheritParams change
#' @param ... Additional arguments passed into method
#'
#' @return A data.frame of numeric values indicating change points and pvalues
#'
#' @noRd
#' @keywords internal
.single_segment <- function(X, method, trim_function, ...) {
  # Trim & stopping criteria
  trim_amt <- trim_function(X)
  nStart <- 1 + trim_amt
  nEnd <- ncol(X) - trim_amt
  if(nStart >= nEnd) return(data.frame('pvalue'=1, 'location'=NA))


  # Find test statistic at every candidate change point
  res <- fchange(X, method=method, type='single', ...)

  # Return location if significant
  data.frame('pvalue'=res$pvalue, 'location'=res$location)
}


#' Change Point Verification
#'
#' This (internal) function is used to verify change points.
#'
#' @inheritParams change
#' @param changes_info A data.frame of numeric values indicating change points and pvalues
#' @param ... Additional arguments passed into method
#'
#' @return A data.frame of numeric values indicating change points and pvalues
#'
#' @noRd
#' @keywords internal
.binary_verification <- function(changes_info, X, method, trim_function,
                                 alpha, silent, ... ) {
  X <- dfts(X)
  if (!silent) cat("-- Verify Step --\n")

  if (!is.null(changes_info)) { # If there was at least one detected
    tmp_changes <- c(0, changes_info$location, ncol(X$data))
    changes_new <- data.frame()
    for (i in 2:(length(tmp_changes) - 1)) {
      ## Get CP
      potential_cp <-
        .single_segment(X=dfts(X=X$data[, (tmp_changes[i - 1] + 1):tmp_changes[i + 1],drop=FALSE],
                               labels = X$labels[(tmp_changes[i - 1] + 1):tmp_changes[i + 1]],
                               fparam = X$fparam, inc.warnings = FALSE),
                        method=method, trim_function=trim_function, ...)

      if(potential_cp$pvalue<=alpha){
        changes_new <- rbind(changes_new,
                         data.frame('location'=potential_cp$location+tmp_changes[i - 1],
                                    'pvalue'=potential_cp$pvalue))
      }
    }
  } else {

    ## Get CP
    changes_new <-
      .single_segment(X=X, method=method, trim_function=trim_function, ...)

    if(changes_new$pvalue>alpha) return()
  }

  # Order and return
  if(nrow(changes_new)<=1) return(changes_new)

  changes_new[order(changes_new$location),]
}
