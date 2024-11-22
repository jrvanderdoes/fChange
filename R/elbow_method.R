#' Elbow Method
#'
#' Method to determine the number of change points using the elbow method. Note,
#'     cascading change points are not considered to allow for every possible
#'     number of change points to be selectable.
#'
#' @param data funts object or Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param test_function Function with the first argument being data. Additional
#'  arguments as passed in as ... . This function should compute a statistic for
#'  each candidate change point (i.e. each curve or column of data), returned as
#'  a named value: allValues.
#' @param trim_function Function with the first element as data (rest with ...).
#'  Used to trim data.
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
#'   },
#'   max_changes=2
#' )
elbow_method <- function(data,
                         test_function = compute_Mn,
                         trim_function = function(...) {
                           0
                         },
                         max_changes = min(ncol(data),20),
                         errorType = "L2",
                         look_ahead = 2, look_alpha = 0.1, ...) {
  data <- .check_data(data)
  max_changes <- min(max_changes,ncol(data))

  # Setup
  n <- ncol(data$data)
  return_data <- data.frame(
    "CP" = NA,
    "Var" = NA
  )

  # Trim & stopping criteria
  trim_amt <- trim_function(data$data, ...)
  if (1 + trim_amt >= n - trim_amt) {
    return()
  }

  # Run First
  stats_tmp <- test_function(data$data[, (1 + trim_amt):(n - trim_amt)], ...)
  test_stat <- c(rep(NA, trim_amt), stats_tmp$allValues, rep(NA, trim_amt))
  tmp_loc <- which.max(stats_tmp$allValues) + trim_amt
  return_data[1, ] <- c(tmp_loc,
                        .compute_total_var(data = data$data, CPs = tmp_loc,
                                           errorType = errorType))

  # Iteratively Search
  while (nrow(return_data) < max_changes) {
    # Setup
    CPs <- c(0, return_data$CP, n)
    CPs <- CPs[order(CPs)]
    prev_CP <- return_data[nrow(return_data), "CP"]
    prev_CP_loc <- which(CPs == prev_CP)

    ## We only need to recompute for the interval changed by last CP!
    # Before
    beforePrevCP <- (CPs[prev_CP_loc - 1] + 1):(CPs[prev_CP_loc])
    trim_amt <- trim_function(data$data[, beforePrevCP], ...)
    test_stat[beforePrevCP] <- NA

    if (1 + trim_amt < length(beforePrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(beforePrevCP) - trim_amt)
      stats_tmp <- test_function(data$data[, beforePrevCP[fill_idx]],
                                 J = length(fill_idx), ...)
      # stats_tmp$allValues[1] <- 0 # Make first zero as don't want to segment
      test_stat[beforePrevCP[fill_idx]] <- stats_tmp$allValues
    }

    # After
    afterPrevCP <- (CPs[prev_CP_loc] + 1):CPs[prev_CP_loc + 1]
    trim_amt <- trim_function(data$data[, afterPrevCP], ...)
    test_stat[afterPrevCP] <- NA

    if (1 + trim_amt < length(afterPrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(afterPrevCP) - trim_amt)
      stats_tmp <- test_function(data$data[, afterPrevCP[fill_idx]],
                                 J = length(fill_idx), ...)
      # stats_tmp$allValues[1] <- 0 # Make first zero as don't want to segment
      test_stat[afterPrevCP[fill_idx]] <- stats_tmp$allValues
    }

    # Temp fix in case of no trim. Will update the check / var calc later
    test_stat <- ifelse(test_stat == 0, NA, test_stat)

    if (sum(is.na(test_stat)) == n) break

    ## Get total variance for each potential CP
    #     TODO: Save previous
    data_segments <- .split_on_NA(test_stat)

    return_data_tmp <- data.frame("CP" = rep(NA,length(data_segments)), "Var" = NA)
    for (k in 1:length(data_segments)) {
      # Find max test stat on interval
      value_max <- max(data_segments[[k]])

      # Get CP and total variance with full data
      section_max <- min(which(test_stat == value_max))
      tmp <- c(section_max, return_data$CP)
      tmp <- tmp[order(tmp)]

      return_data_tmp[k, ] <- c(
        section_max,
        .compute_total_var(data$data, tmp, errorType)
      )
    }

    ## With max test-statistic on each section, take one leading to min variance
    return_data[nrow(return_data)+1, ] <- return_data_tmp[which.min(return_data_tmp$Var), ]
  }

  # Add No Change option and compute percent change
  return_data <- rbind(
    c(NA, .compute_total_var(data$data, c(), errorType)),
    return_data
  )
  return_data$Percent <- 1 - return_data$Var / max(return_data$Var)
  return_data$Gain <- c(NA,return_data$Var[-1] /return_data$Var[-nrow(return_data)])

  # Define vars to remove notes
  CP <- Var <- Percent <- NULL

  var_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Total Variance") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  per_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent),size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Percent Explained") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  gain_plot <-
    ggplot2::ggplot(return_data[-1,]) +
    ggplot2::geom_point(ggplot2::aes(x = 1:length(CP), y = Gain), size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 1:length(CP), y = Gain),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Percent Explained") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  data_change<- data.frame(matrix(NA,nrow=nrow(return_data),ncol = look_ahead))
  for(i in 1:look_ahead){
    data_change[,i] <- c(return_data$Gain[-c(1:i)]>1-look_alpha,rep(TRUE,i))
  }
  cutoff <- which.max(apply(data_change, MARGIN = 1, prod)) - 1

  recommend_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),
                       linewidth=2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=cutoff),
                        linetype='dotted', col='red',linewidth=2) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),
                        size=4) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Total Variance") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  list("CPInfo" = return_data, "VariancePlot" = var_plot,
       "PercentPlot" = per_plot, "GainPlot" = gain_plot,
       "RecommendPlot" = recommend_plot)
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
  CPs <- unique(c(0, CPs, ncol(data)))

  if (errorType == "L2" || errorType == "Tr") {
    # Standardize Data
    for (i in 2:length(CPs)) {
      indices <- (CPs[i - 1] + 1):CPs[i]
      CP_mean <- rowMeans(as.data.frame(data[, indices])) ## TODO:: Verify
      # Standardize Data
      data_std[, indices] <- data_std[, indices] - CP_mean
    }

    # ## Cannot use cov because it removes the mean from individual FDs, thus
    # #     making each change point near equally effective
    # covMatrix <- matrix(NA, ncol = ncol(data), nrow = ncol(data))
    # sampleCoef <- 1 / (nrow(data) - 1)
    #
    # for (i in 1:ncol(data)) {
    #   for (j in i:ncol(data)) {
    #     covMatrix[i, j] <- covMatrix[j, i] <-
    #       sum(sampleCoef * (data_std[, i] * data_std[, j]))
    #   }
    # }
    #
    covMatrix <- autocov_approx_h(data_std,0)

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
