#' Title
#'
#' @param data_length
#' @param dep
#' @param nSims
#' @param path
#'
#' @return
#'
#' @examples
#' simulate_data_size(100,0)
#' simulate_data_size(100,0.5)
simulate_data_size <- function(data_length = 100, dep = 0, nSims = 500,
                               path = tempdir()) {
  data <- list()
  for (i in 1:nSims) {
    cat(paste0(i, ", "))
    set.seed(data_length * 10000 + dep * 1000 + i)
    data[[i]] <- generate_data_fd(
      ns = c(data_length),
      eigsList = list(c(3, 2, 1, 0.5)),
      basesList = list(fda::create.bspline.basis(nbasis = 4, norder = 4)),
      meansList = c(0),
      distsArray = c("Normal"),
      evals = seq(0, 1, 0.05),
      kappasArray = c(dep),
      silent = T
    )
  }


  saveRDS(data, paste0(path, "data", "_dl", data_length, "_d", dep, ".rds"))
}


#' Title
#'
#' @param data_length
#' @param dep
#' @param nSims
#' @param path
#'
#' @return
#'
#' @examples
#' simulate_size(100,0)
#' simulate_size(100,0.5)
simulate_size <- function(data_length = 100, dep = 0, nSims = 500,
                          path = tempdir()) {
  result <- data.frame(
    "Cov_iid" = rep(NA, nSims),
    "Cov_fn3" = NA,
    "Cov_f2n3" = NA,
    "Cov_f3n3" = NA,
    "Perm_iid" = NA,
    "Perm_fn3" = NA,
    "Perm_f2n3" = NA,
    "Perm_f3n3" = NA,
    "Welch_iid" = NA,
    "Welch_fn3" = NA,
    "Welch_f2n3" = NA,
    "Welch_f3n3" = NA
  )
  data <- readRDS(paste0(path, "data", "_dl", data_length, "_d", dep, ".rds"))
  for (i in 1:nSims) {
    cat(paste0(i, ", "))
    set.seed(data_length * 10000 + dep * 1000 + i + 1) # +1 to change from data seeds
    Cov_iid <- detect_changepoint_singleCov(
      X = data[[i]], nSims = 1000, h = 0,
      space = "BM", Cov_M = 20, silent = T
    )
    Cov_n3 <- detect_changepoint_singleCov(
      X = data[[i]], nSims = 1000, h = data_length^(1 / 3),
      space = "BM", Cov_M = 20, silent = T
    )
    Cov_2n3 <- detect_changepoint_singleCov(
      X = data[[i]], nSims = 1000, h = 2 * data_length^(1 / 3),
      space = "BM", Cov_M = 20, silent = T
    )
    Cov_3n3 <- detect_changepoint_singleCov(
      X = data[[i]], nSims = 1000, h = 3 * data_length^(1 / 3),
      space = "BM", Cov_M = 20, silent = T
    )

    Perm_iid <- generalized_resampling(
      X = data[[i]], blockSize = 1, fn = compute_Tn,
      space = "BM", M = 1000, silent = T
    )
    Perm_n3 <- generalized_resampling(
      X = data[[i]], blockSize = data_length^(1 / 3),
      fn = compute_Tn, space = "BM", M = 1000, silent = T
    )
    Perm_2n3 <- generalized_resampling(
      X = data[[i]], blockSize = 2 * data_length^(1 / 3),
      fn = compute_Tn, space = "BM", M = 1000, silent = T
    )
    Perm_3n3 <- generalized_resampling(
      X = data[[i]], blockSize = 3 * data_length^(1 / 3),
      fn = compute_Tn, space = "BM", M = 1000, silent = T
    )

    Approx_iid <- welch_approximation(
      X = data[[i]], alpha = 0.05, TVal = ncol(data[[i]]),
      W = NULL, W1 = NULL, M = 1000, h = 0,
      K = bartlett_kernel
    )
    Approx_n3 <- welch_approximation(
      X = data[[i]], alpha = 0.05, TVal = ncol(data[[i]]),
      W = NULL, W1 = NULL, M = 1000, h = data_length^(1 / 3),
      K = bartlett_kernel
    )
    Approx_2n3 <- welch_approximation(
      X = data[[i]], alpha = 0.05, TVal = ncol(data[[i]]),
      W = NULL, W1 = NULL, M = 1000, h = 2 * data_length^(1 / 3),
      K = bartlett_kernel
    )
    Approx_3n3 <- welch_approximation(
      X = data[[i]], alpha = 0.05, TVal = ncol(data[[i]]),
      W = NULL, W1 = NULL, M = 1000, h = 3 * data_length^(1 / 3),
      K = bartlett_kernel
    )

    Tn <- mean(Cov_iid$value, Cov_n3$value, Cov_3n3$value)

    result[i, ] <- c(
      Cov_iid$pval <= 0.05,
      Cov_n3$pval <= 0.05,
      Cov_2n3$pval <= 0.05,
      Cov_3n3$pval <= 0.05,
      Perm_iid$pval <= 0.05,
      Perm_n3$pval <= 0.05,
      Perm_2n3$pval <= 0.05,
      Perm_3n3$pval <= 0.05,
      Approx_iid <= Tn,
      Approx_n3 <= Tn,
      Approx_2n3 <= Tn,
      Approx_3n3 <= Tn
    )
    saveRDS(
      result,
      paste0(
        path, "results",
        "_dl", data_length, "_d", dep, ".rds"
      )
    )
  }

  result
}

