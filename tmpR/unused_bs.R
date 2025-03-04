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
#' @examples
#' #.single_binary_segmentation(
#' #   data = electricity[, 1:60],
#' #   test_statistic_function = compute_Mn,
#' #   cutoff_function = welch_approximation,
#' #   trim_function = function(data) {
#' #     max(10, floor(log(ncol(as.data.frame(data)))),
#' #     na.rm = TRUE
#' #  )
#' # }
#' #)
.single_binary_segmentation <- function(data,
                                        test_statistic_function,
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
  test_stats <- test_statistic_function(as.data.frame(data))

  # Return index of max change point if larger than cutoff
  return_value <- ifelse(test_stats$value >= cutoff_function(data, alpha=alpha, ...),
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
#'     \code{.single_binary_segmentation}
#'
#' @return Numeric values indicating the change points detected
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup Data
#' data_KL <- generate_kl(
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
#'   data_KL$data, compute_Tn, welch_approximation,
#'   function(data) {
#'     max(2, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   }
#' )
#' generalized_wild_binary_segmentation(
#'   data = data_KL$data,
#'   test_statistic_function = compute_Tn,
#'   cutoff_function = welch_approximation,
#'   trim_function = function(data) {
#'     max(2, floor(log(ncol(as.data.frame(data)))),
#'       na.rm = TRUE
#'     )
#'   }
#' )
#' }
generalized_wild_binary_segmentation <- function(data, M = 5000, add_full = TRUE, block_size = 1,
                                                 ...) {
  # Setup
  n <- ncol(data)
  changes <- c()
  result <- matrix(ncol = 2, nrow = M + add_full)

  # Test
  if (n <= 1) {
    return()
  }

  # Run
  if (add_full) {
    result[1, ] <- .single_binary_segmentation(data, include_value = TRUE, ...)
  }

  for (i in add_full + 1:M) {
    min_pt <- sample(1:(n - block_size + 1), 1)
    max_pt <- sample((min_pt + block_size - 1):n, 1)

    # Must return location and value
    # This needs to return location and values!
    result[i, ] <- .single_binary_segmentation(data[, min_pt:max_pt],
                                               include_value = TRUE, ...
    )
  }

  # Select best and continue if reasonable
  if (nrow(stats::na.omit(result))) {
    cp_loc <- result[which.max(result[, 2]), 1]

    changes <- c(
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

  changes
}
