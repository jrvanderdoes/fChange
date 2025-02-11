#' Complete Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. Change points are recursively found until no
#'     more change points are detected.
#'
#' @param X funts object or numeric data.frame with rows for evaluated values and columns
#'    indicating functional observations
#' @param test_statistic_function XXXXXXXXXXXXXXXXX.
#' @param cutoff_function XXXXXX
#' @param trim_function XXXXXX
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
#' @param ... Additional arguments passed into test_statistic_function
#'
#' @return A list of numeric values indicating change points  (if exists),
#'     NA otherwise
#'
#' @noRd
#' @keywords internal
.binary_segmentation <- function(X, method,
                                 trim_function = function(X) {
                                   max(10, floor(log(ncol(as.data.frame(X)))),
                                       na.rm = TRUE
                                   )
                                 },
                                 alpha = 0.05, silent = FALSE, ...) {
  ## Setup
  X <- funts(X)

  ## Get change points
  CPsVals <- .detect_changes(
    X = X, method = method,
    trim_function = trim_function,
    alpha = alpha,
    addAmt = 0,
    silent = silent,
    ...
  )

  ## Verify
  CPsVals <- .binary_verification(
    CPsVals = CPsVals, X = X, method = method,
    trim_function = trim_function, alpha = alpha,
    silent = silent, ...
  )

  ## Return Results
  CPsVals
}



#' Detect Change Points
#'
#' This (internal) function is multiple .single_binary_segmentation for
#'     complete_binary_segmentation. It recursively calls itself
#'
#' @param X Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#'
#' @param test_statistic_function Function with the first argument being data.
#'     Additional arguments passed in via ... . Return the selected change point.
#'     Give this or test_statistic_function.
#' @param alpha (Optional) Numeric value in \[0,1\] indicating the significance.
#' @param addAmt (Optional) Default 0. This adds a number to shift change points.
#'               Used for recursive calls, likely no need to change.
#' @param silent (Optional) Boolean indicating if output should be given. Default
#'               is FALSE (meaning print output).
#' @param ... (Optional) Additional parameters for the functions.
#'
#' @return Vector of detected change point locations
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
  if (potential$pvalue > alpha) {
    return()
  }

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
      X = funts(X$data[, 1:potential$location,drop=FALSE],
                labels = X$labels[1:potential$location],
                intraobs = X$intraobs,inc.warnings = FALSE),
      method = method,
      trim_function = trim_function,
      alpha = alpha,
      addAmt = addAmt,
      silent = silent,
      ...
    ),
    data.frame('location'=potential$location + addAmt, 'pvalue'=potential$pvalue),
    .detect_changes(
      X = funts(X$data[, (potential$location + 1):ncol(X)],labels = X$labels[(potential$location + 1):ncol(X)],
                        intraobs = X$intraobs,inc.warnings = FALSE),
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
#' @inheritParams generalized_binary_segmentation
#' @param trim_function XXXXX
#' @param include_value XXXXX
#'
#' @return A numeric value indicating the cutoff location (if exists),
#'     NA otherwise
#'
#' @noRd
#' @keywords internal
.single_segment <- function(X, method, trim_function, ...) {
  # Trim & stopping criteria
  trim_amt <- trim_function(X)
  nStart <- 1 + trim_amt
  nEnd <- ncol(X) - trim_amt
  if(nStart >= nEnd) return(NA)


  # Find test statistic at every candidate change point
  res <- change(X, method=method, type='single', ...)

  # Return location if significant
  data.frame('pvalue'=res$pvalue, 'location'=res$location)
}


#' Change Point Verification
#'
#' This (internal) function is used to verify change points.
#'
#' @param data funts object or numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param CPsVals Numeric vector indicating change point locations (empty vector
#'     used if no change point detected)
#' @param test_statistic_function Function with the first argument being data
#'     and the second argument optional argument for candidate change points.
#'     Additional arguments passed in via ... . Return a single numeric value.
#' @param cutoff_function Function with first argument being data and the second
#'     argument being alpha. No other arguments given. Return single numeric
#'     value.
#' @param trim_function Function taking data as an argument and returning a
#'     numeric value indicating how much should be trimmed on each end
#' @param alpha Numeric value in \[0,1\] indicating the significance for
#'     cutoff_function.
#' @param silent Boolean to indicate if progress output should be printed
#' @param ... Additional inputs to pass to the given functions
#'
#' @return CPsVals Numeric vector indicating change point locations (NA if no
#'     change points are detected)
#'
#' @noRd
#' @keywords internal
.binary_verification <- function(CPsVals, X, method, trim_function,
                                 alpha, silent, ... ) {
  X <- funts(X)
  if (!silent) cat("-- Verify Step --\n")

  if (!is.null(CPsVals)) { # If there was at least one detected
    tmp_cps <- c(0, CPsVals$location, ncol(X$data))
    CPsVals_new <- data.frame()
    for (i in 2:(length(tmp_cps) - 1)) {
      ## Get CP
      potential_cp <-
        .single_segment(X=funts(X=X$data[, (tmp_cps[i - 1] + 1):tmp_cps[i + 1]],
                                labels = X$labels[(tmp_cps[i - 1] + 1):tmp_cps[i + 1]],
                                intraobs = X$intraobs, inc.warnings = FALSE),
                        method=method, trim_function=trim_function, ...)
      if(potential_cp$pvalue<=alpha){
        CPsVals_new <- rbind(CPsVals_new,
                         data.frame('location'=potential_cp$location+tmp_cps[i - 1],
                                    'pvalue'=potential_cp$pvalue))
      }
    }
  } else {

    ## Get CP
    CPsVals_new <-
      .single_segment(X=X, method=method, trim_function=trim_function, ...)

    if(CPsVals_new$pvalue>alpha) return()
  }

  # Order and return
  if(nrow(CPsVals_new)<=1) return(CPsVals_new)

  CPsVals_new[order(CPsVals_new$location),]
}
