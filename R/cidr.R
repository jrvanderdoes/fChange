#' Compute (Overnight) Cumulative Intraday Returns
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#'
#' @returns A dfts object with CIDR or OCIDRs
#' @export
#'
#' @references Rice, G., Wirjanto, T., & Zhao, Y. (2023). Exploring volatility
#'  of crude oil intraday return curves: A functional GARCH-X model. Journal of
#'  Commodity Markets, 32, 100361-.
#'
#' @examples
#' tmp <- dfts(SPYUS500$data[, 1:100],
#'   name = "SP500 100 Days",
#'   labels = SPYUS500$labels[1:100], fparam = SPYUS500$fparam
#' )
#' cidr(tmp)
cidr <- function(X) {
  X <- dfts(X)

  dat_cidr <- X$data
  for (i in 1:nrow(dat_cidr)) {
    dat_cidr[i, ] <- 100 * (log(X$data[i, ]) - log(X$data[1, ]))
  }

  dfts(dat_cidr,
    name = paste0("CIDR of ", X$name), labels = X$labels,
    fparam = X$fparam, inc.warnings = F
  )
}


#' @rdname cidr
#'
#' @export
#'
#' @examples
#' tmp <- dfts(SPYUS500$data[, 1:100],
#'   name = "SP500 100 Days",
#'   labels = SPYUS500$labels[1:100], fparam = SPYUS500$fparam
#' )
#' ocidr(tmp)
ocidr <- function(X) {
  X <- dfts(X)

  dat_cidr <- X$data
  dat_cidr[, 1] <- 100 * (log(dat_cidr[, 1]) - log(dat_cidr[1, 1]))
  # Obs
  for (j in 2:ncol(dat_cidr)) {
    dat_cidr[, j] <- 100 *
      (log(X$data[, j]) - log(X$data[nrow(dat_cidr), j - 1]))
  }

  dfts(dat_cidr,
    name = paste0("OCIDR of ", X$name), labels = X$labels,
    fparam = X$fparam, inc.warnings = F
  )
}
