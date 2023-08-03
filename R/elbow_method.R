#' Elbow Method
#'
#' Method to determine the number of change points using the elbow method. Note,
#'     cascading change points are not considered to allow for every possible
#'     number of change points to be selectable.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param test_function Function with the first argument being data. Additional
#'  arguments as passed in as ... . This function should compute a statistic for
#'  each candidate change point (i.e. each curve or column of data), returned as
#'  a named value: allValues.
#' @param trim_function XXXXXXX
#' @param errorType String of 'L2' or 'Tr' indicating the error function to use
#' @param ... Additional parameters to pass into the respective functions
#'
#' @return list with three elements.
#'  \itemize{
#'    \item **CPInfo**: Data.frame of candidate changes, ordered by impact. The
#'      columns are: CP, Var, and Percent. CP indicates change location, Var is
#'      total variance, Percent is the percentage of variance removed compared
#'      to the previous step. The first row has CP=NA, indicating the no change
#'      scenario.
#'    \item **VarPlot**: A ggplot object showing the variance for increasing
#'      number of changes.
#'    \item **PerPlot**: A ggplot object showing the percent of variance removed
#'      compared to the previous step.
#'  }
#' @export
#'
#' @examples
#' results <- elbow_method(
#'   data = electricity[, 1:50],
#'   trim_function = function(...) {
#'     10
#'   }
#' )
elbow_method <- function(data,
                         test_function = compute_Mn,
                         trim_function = function(...) {
                           0
                         },
                         errorType = "L2", ...) {
  # Setup
  n <- ncol(data)
  return_data <- data.frame(
    "CP" = NA,
    "Var" = NA
  )

  # Trim & stopping criteria
  trim_amt <- trim_function(data, ...)
  if (1 + trim_amt > n - trim_amt) {
    return()
  }

  # Run First
  stats_tmp <- test_function(data[, (1 + trim_amt):(n - trim_amt)], ...)
  test_stat <- c(rep(NA, trim_amt), stats_tmp$allValues, rep(NA, trim_amt))
  tmp_loc <- which.max(stats_tmp$allValues) + trim_amt
  return_data[1, ] <- c(tmp_loc, .compute_total_var(data, tmp_loc, errorType))

  # Iteratively Search
  i <- 1
  while (TRUE) {
    # Setup
    i <- i + 1
    CPs <- c(0, return_data$CP, n)
    CPs <- CPs[order(CPs)]
    prev_CP <- return_data[nrow(return_data), "CP"]
    prev_CP_loc <- which(CPs == prev_CP)

    ## We only need to recompute for the interval changed by last CP!
    # Before
    beforePrevCP <- (CPs[prev_CP_loc - 1] + 1):(CPs[prev_CP_loc])
    trim_amt <- trim_function(data[, beforePrevCP], ...)
    test_stat[beforePrevCP] <- NA

    if (1 + trim_amt < length(beforePrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(beforePrevCP) - trim_amt)
      stats_tmp <- test_function(data[, beforePrevCP[fill_idx]], ...)
      test_stat[beforePrevCP[fill_idx]] <- stats_tmp$allValues
    }

    # After
    afterPrevCP <- (CPs[prev_CP_loc] + 1):CPs[prev_CP_loc + 1]
    trim_amt <- trim_function(data[, afterPrevCP], ...)
    test_stat[afterPrevCP] <- NA

    if (1 + trim_amt < length(afterPrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(afterPrevCP) - trim_amt)
      stats_tmp <- test_function(data[, afterPrevCP[fill_idx]], ...)
      test_stat[afterPrevCP[fill_idx]] <- stats_tmp$allValues
    }

    # Temp fix in case of no trim. Will update the check / var calc later
    test_stat <- ifelse(test_stat == 0, NA, test_stat)

    if (sum(is.na(test_stat)) == n) break

    ## Get total variance for each potential CP
    #     TODO: Save previous
    data_segments <- .split_on_NA(test_stat)

    return_data_tmp <- data.frame("CP" = NA, "Var" = NA)
    for (k in 1:length(data_segments)) {
      # Find max test stat on interval
      value_max <- max(data_segments[[k]])

      # Get CP and total variance with full data
      section_max <- which(test_stat == value_max)
      tmp <- c(section_max, return_data$CP)
      tmp <- tmp[order(tmp)]

      return_data_tmp[k, ] <- c(
        section_max,
        .compute_total_var(data, tmp, errorType)
      )
    }

    ## With max test-statistic on each section, take one leading to min variance
    return_data[i, ] <- return_data_tmp[which.min(return_data_tmp$Var), ]


    if (nrow(return_data) == (n - 1)) break
  }

  # Add No Change option and compute percent change
  return_data <- rbind(
    c(NA, .compute_total_var(data, c(), errorType)),
    return_data
  )
  return_data$Percent <- 1 - return_data$Var / max(return_data$Var)

  # Define vars to remove notes
  CP <- Var <- Percent <- NULL

  var_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Var)) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Var)) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Total Variance") +
    ggplot2::theme_bw()

  per_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent)) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent)) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Percent Explained") +
    ggplot2::theme_bw()

  list("CPInfo" = return_data, "VarPlot" = var_plot, "PerPlot" = per_plot)
}

#' Compute Total Variance
#'
#' This (internal) function computes the total variance in the data with given
#'  CPs.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD.
#' @param CPs Vector of numerics indicating the changepoint locations
#' @param errorType (Optional) String of 'L2' or 'Tr' indicating the error
#'                  function to use. Default is L2.
#' @param M (Optional) Numeric indicating the number of Brownian motions needed
#'          for CE error.
#'
#' @return Numeric indicating the variance between all subsegments
#'
#' @noRd
.compute_total_var <- function(data, CPs, errorType = "L2", M = 1000) {
  # Setup
  data <- as.data.frame(data)
  data_std <- data
  CPs <- c(0, CPs, ncol(data))

  if (errorType == "L2" || errorType == "Tr") {
    # Standardize Data
    for (i in 2:length(CPs)) {
      indices <- (CPs[i - 1] + 1):CPs[i]
      CP_mean <- rowMeans(as.data.frame(data[, indices])) ## TODO:: Verify
      # Standardize Data
      data_std[, indices] <- data_std[, indices] - CP_mean
    }

    ## Cannot use cov because it removes the mean from individual FDs, thus
    #     making each change point near equally effective
    covMatrix <- matrix(NA, ncol = ncol(data), nrow = ncol(data))
    sampleCoef <- 1 / (nrow(data) - 1)

    for (i in 1:ncol(data)) {
      for (j in i:ncol(data)) {
        covMatrix[i, j] <- covMatrix[j, i] <-
          sum(sampleCoef * (data_std[, i] * data_std[, j]))
      }
    }
  } else if (errorType == "CE") {
    W <- as.data.frame(sapply(rep(0, M), sde::BM, N = nrow(data) - 1))
    CE <- data.frame(matrix(ncol = ncol(data), nrow = M))

    ## Compute CEs
    for (i in 1:ncol(data)) {
      CE[, i] <- apply(W, 2,
        FUN = function(v, dat) {
          exp(complex(real = 0, imaginary = 1) * (t(dat) %*% v))
        }, dat = data[, i]
      )
    }

    CE_std <- CE

    # Standardize Data
    for (i in 2:length(CPs)) {
      indices <- (CPs[i - 1] + 1):CPs[i]
      CP_mean <- rowMeans(as.data.frame(CE[, indices])) ## TODO:: Verify
      # Standardize Data
      CE_std[, indices] <- CE[, indices] - CP_mean
    }

    # covMatrix <- matrix(NA,ncol=ncol(data), nrow = ncol(data))
    # sampleCoef <- 1/(nrow(data)-1)
    #
    # ## Compute covMatrix
    # for(i in 1:ncol(data)){
    #   for(j in (i):ncol(data)){
    #     covMatrix[i,j] <- covMatrix[j,i] <-
    #       sum(sampleCoef*(CE[,i] %*% CE[,j]))
    #   }
    # }
  }

  # Apply error metric
  if (errorType == "L2") {
    returnValue <- sqrt(sum(covMatrix^2))
  } else if (errorType == "Tr") {
    returnValue <- sum(diag(covMatrix))
  } else if (errorType == "CE") {
    returnValue <- sum(colSums(abs(CE_std)^2) / M) # sum(colSums(abs(CE_std)^2)/M)
  } else {
    stop("Only L2, Tr, and CE error functions implemented")
  }

  returnValue
}


#' Split on NA
#'
#' This (internal) function splits a vector based on any NA values.
#'
#' @param vec Vector to be split containing numerics and NA values.
#'
#' @return List with each item being a group separated by NAs.
#'
#' @noRd
.split_on_NA <- function(vec) {
  is.sep <- is.na(vec)
  split(vec[!is.sep], cumsum(is.sep)[!is.sep])
}
