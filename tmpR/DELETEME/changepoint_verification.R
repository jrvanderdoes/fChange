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
.changepoint_verification <- function(data, CPsVals,
                                      test_statistic_function = NULL,
                                      changepoint_function = NULL,
                                      silent = FALSE,
                                      ...) {
  data <- funts(data)
  if (!silent) cat("-- Verify Step --\n")

  if (length(CPsVals) >= 1) { # If there was a CP detected
    tmp_cps <- c(0, CPsVals, ncol(data$data))
    tmp_cps <- tmp_cps[order(tmp_cps)]
    CPsVals <- c()
    for (i in 2:(length(tmp_cps) - 1)) {
      # Get CP
      if (!is.null(test_statistic_function)) {
        potential_cp <-
          .single_binary_segmentation(data$data[, (tmp_cps[i - 1] + 1):tmp_cps[i + 1]],
            test_statistic_function = test_statistic_function, ...
          )
      } else if (!is.null(changepoint_function)) {
        potential_cp <- changepoint_function(data$data[, (tmp_cps[i - 1] + 1):tmp_cps[i + 1]], ...)
      }

      if (!is.na(potential_cp)) {
        CPsVals <- c(CPsVals, potential_cp + tmp_cps[i - 1])
      }
    }
  } else {
    # Get CP
    # TODO:: Check this out
    if (!is.null(test_statistic_function)) {
      CPsVals <-
        .single_binary_segmentation(data$data,
          test_statistic_function = test_statistic_function, ...
        )
    } else if (!is.null(changepoint_function)) {
      CPsVals <- changepoint_function(data$data, ...)
    }
  }

  # Order and return
  if (sum(is.na(CPsVals)) == length(CPsVals)) {
    return(NA)
  }
  CPsVals[order(CPsVals)]
}
