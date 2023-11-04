#' Complete Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. Change points are recursively found until no
#'     more change points are detected.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating functional observations
#' @param changepoint_function XXXXXXXXXXXXXXXXX.
#' @param cutoff_function XXXXXX
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
#' @param ... Additional arguments passed into changepoint_function
#'
#' @return A list of numeric values indicating change points  (if exists),
#'     NA otherwise
#' @export
#'
#' @examples
#' # Below will give NA as the trim cuts too much.
#' #   change 50 to 20 in trim function to get other results
#' complete_binary_segmentation(
#'   data = electricity[, 1:80],
#'   changepoint_function = compute_Mn,
#'   cutoff_function = welch_approximation,
#'   trim_function = function(data) {
#'     max(50, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   },
#'   final_verify = FALSE
#' )
complete_binary_segmentation <- function(data,
                                         changepoint_function = compute_Mn,
                                         cutoff_function = welch_approximation,
                                         final_verify = TRUE,
                                         silent = FALSE,
                                         alpha = 0.05,
                                         ...) {
  # Get change points
  CPsVals <- .detectChangePoints(
    data = data,
    changepoint_function = changepoint_function,
    cutoff_function = cutoff_function,
    silent = silent, alpha = alpha, ...
  )

  # Verify as desired
  if (final_verify) {
    CPsVals <- .changepoint_verification(
      CPsVals = CPsVals, data = data, changepoint_function = changepoint_function,
      cutoff_function = cutoff_function, silent = silent, alpha = alpha, ...
    )
  }

  return(CPsVals)
}


#' Single Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. At most one change point is detected.
#'
#' @inheritParams complete_binary_segmentation
#' @param trim_function XXXXX
#' @param include_value XXXXX
#'
#' @return A numeric value indicating the cutoff location (if exists),
#'     NA otherwise
#' @export
#'
#' @examples
#' single_binary_segmentation(
#'   data = electricity[, 1:60],
#'   changepoint_function = compute_Mn,
#'   cutoff_function = welch_approximation,
#'   trim_function = function(data) {
#'     max(10, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   }
#' )
single_binary_segmentation <- function(data, changepoint_function,
                                       cutoff_function,
                                       trim_function,
                                       alpha = 0.05, include_value = FALSE,
                                       ...) {
  # Trim & stopping criteria
  trim_amt <- trim_function(data, ...)
  nStart <- 1 + trim_amt
  nEnd <- ncol(as.data.frame(data)) - trim_amt
  if (nStart >= nEnd) ifelse(include_value, return(c(NA, NA)), return(NA))

  # Find test statistic at every candidate change point
  test_stats <- changepoint_function(as.data.frame(data))

  # Return index of max change point if larger than cutoff
  return_value <- ifelse(test_stats$value >= cutoff_function(data, alpha, ...),
    test_stats$location,
    NA
  )

  # Add in value
  if (include_value) {
    if (!is.na(return_value)) {
      return_value <- c(return_value, max(test_stats, na.rm = TRUE))
    } else {
      return_value <- c(return_value, NA)
    }
  }

  return_value
}


#' Wild Binary Segmentation
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param M Numeric value
#'
#' Indicates the number of intervals to examine
#'
#' @param add_full Boolean value
#'
#' Indicates if the entire interval should be added
#'
#' @param block_size (Optional) Numeric value
#'
#' Indicates the minimum block size
#'
#' @param ... Additional parameters to pass into
#'     \code{single_binary_segmentation}
#'
#' @return Numeric values indicating the change points detected
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup Data
#' data_KL <- generate_data_fd(
#'   ns = c(12, 12, 12),
#'   eigsList = list(
#'     c(3, 2, 1, 0.5),
#'     c(3, 2, 1, 0.5),
#'     c(3, 2, 1, 0.5)
#'   ),
#'   basesList = list(
#'     fda::create.bspline.basis(nbasis = 4, norder = 4),
#'     fda::create.bspline.basis(nbasis = 4, norder = 4),
#'     fda::create.bspline.basis(nbasis = 4, norder = 4)
#'   ),
#'   meansList = c(-1, 0, 1),
#'   distsArray = c("Normal", "Normal", "Normal"),
#'   evals = seq(0, 1, 0.05),
#'   kappasArray = c(0, 0, 0)
#' )
#'
#' complete_binary_segmentation(
#'   data_KL, compute_Tn, welch_approximation,
#'   function(data) {
#'     max(2, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   }
#' )
#' wild_binary_segmentation(
#'   data = data_KL,
#'   test_statistic_function = compute_Tn,
#'   cutoff_function = welch_approximation,
#'   trim_function = function(data) {
#'     max(2, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   }
#' )
#' }
wild_binary_segmentation <- function(data, M = 5000, add_full = TRUE, block_size = 1,
                                     ...) {
  # Setup
  n <- ncol(data)
  cps <- c()
  result <- matrix(ncol = 2, nrow = M + add_full)

  # Test
  if (n <= 1) {
    return()
  }

  # Run
  if (add_full) {
    result[1, ] <- single_binary_segmentation(data, include_value = TRUE, ...)
  }

  for (i in add_full + 1:M) {
    min_pt <- sample(1:(n - block_size + 1), 1)
    max_pt <- sample((min_pt + block_size - 1):n, 1)

    # Must return location and value
    # This needs to return location and values!
    result[i, ] <- single_binary_segmentation(data[, min_pt:max_pt],
      include_value = TRUE, ...
    )
  }

  # Select best and continue if reasonable
  if (nrow(stats::na.omit(result))) {
    cp_loc <- result[which.max(result[, 2]), 1]

    cps <- c(
      wild_binary_segmentation(data[, min_pt:cp_loc],
        M = M, add_full = add_full, block_size = block_size,
        ...
      ),
      cp_loc,
      wild_binary_segmentation(data[, (cp_loc + 1):max_pt],
        M = M, add_full = add_full, block_size = block_size,
        ...
      ) + cp_loc
    )
  }

  cps
}


#' Detect Change Points
#'
#' This (internal) function is multiple single_binary_segmentation for
#'     complete_binary_segmentation. It recursively calls itself
#'
#' TODO: Remove alpha
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#'
#' @param changepoint_function Function with the first argument being data.
#'     Additional arguments passed in via ... . Return the selected change point.
#'     Give this or changepoint_function.
#' @param alpha (Optional) Numeric value in [0,1] indicating the significance.
#' @param addAmt (Optional) Default 0. This adds a number to shift change points.
#'               Used for recursive calls, likely no need to change.
#' @param silent (Optional) Boolean indicating if output should be given. Default
#'               is FALSE (meaning print output).
#' @param ... (Optional) Additional parameters for the functions.
#'
#' @return Vector of detected change point locations
#'
#' @noRd
.detectChangePoints <- function(data,
                                changepoint_function,
                                alpha = NULL,
                                addAmt = 0,
                                silent = FALSE,
                                ...) {
  # Look for a single change
  potential_cp <- single_binary_segmentation(data,
    changepoint_function = changepoint_function, alpha = alpha, ...
  )

  # No Change Point Detected
  if (is.na(potential_cp)) {
    return()
  }

  # Display progress
  if (!silent) {
    cat(paste0(
      "ChangePoint Detected (", 1 + addAmt, "-", addAmt + ncol(data), " at ",
      addAmt + potential_cp, "): Segment Data and Re-Search\n"
    ))
  }

  # Search Recursively
  return(c(
    .detectChangePoints(
      data = as.data.frame(data[, 1:potential_cp]),
      changepoint_function = changepoint_function,
      addAmt = addAmt,
      silent = silent,
      alpha = alpha,
      ...
    ),
    potential_cp + addAmt,
    .detectChangePoints(
      data = as.data.frame(data[, (potential_cp + 1):ncol(data)]),
      changepoint_function = changepoint_function,
      addAmt = addAmt + potential_cp,
      silent = silent,
      alpha = alpha,
      ...
    )
  ))
}
