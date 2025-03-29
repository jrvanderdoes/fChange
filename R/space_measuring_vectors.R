#' Compute Spacing Measuring Functions
#'
#' This function is used to compute discretized functions, i.e. vectors, to
#'  explore functional spaces.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param M Integer for the number of functions to generate.
#' @param space String for the space of interest. Options are Brownian motion
#' ('BM'), principal components ('PC'), and vectors in iid standard,
#' random normals ('RN'). Additional options are forthcoming
#'
#' @return Data.frame with columns of discretized functions describing the space.
#'  Columns are independent functions.
#' @export
#'
#' @seealso [fchange()]
#'
#' @examples
#' space_measuring_functions(M=10, space="BM", X=electricity)
#' space_measuring_functions(M=10, space="PC", X=electricity)
space_measuring_functions <- function(X, M=20, space='BM') {
  X <- dfts(X)

  if (space == "BM") {
    # W <- as.data.frame(sapply(rep(0, M), sde::BM, N = nrow(X$data) - 1))
    W <- generate_brownian_motion(M,v=seq(0,1,length.out=nrow(X$data)))$data
  } else if (space == "PC") {
    pComps <- stats::prcomp(X$data, center = FALSE, scale = FALSE)
    W <- as.data.frame(sapply(rep(nrow(X$data), M),
      function(x, pcs) {
        rowSums(stats::rnorm(ncol(pcs)) * pcs)
      },
      pcs = pComps$x
    ))
  } else if (space == "OU") {
    # TODO:: Fix this
    stop('Sorry, OU space is under testing.',call. = FALSE)
    x <- seq(0, 1, length.out = nrow(X$data))
    covMat <- (matrix(1, ncol = length(x), nrow = length(x)))
    covMat[, 1] <- covMat[1, ] <- exp(abs(x - x[1]))
    for (i in 2:(length(x) - 1)) {
      covMat[, i] <- covMat[i, ] <- c(
        covMat[i, c(1:(i - 1))],
        exp(abs(x - x[i]))[-c(1:(i - 1))]
      )
    }
    W <- as.data.frame(sapply(rep(nrow(X$data), M),
      function(x, covMat) {
        covMat %*% stats::rnorm(x)
      },
      covMat = covMat
    ))

    # W <- as.data.frame(sapply(rep(nrow(X),M),sde::rsOU,theta=c(0,0,1)))
  } else if (space == "RN") {
    W <- data.frame(matrix(1, ncol = M, nrow = nrow(X$data) - 1))
    W <- as.data.frame(sapply(
      rep(nrow(X$data), M),
      function(x) {
        rep(stats::rnorm(1), x)
      }
    ))
  } else {
    stop("Error: Sorry only BM, OU, PC, or RN processes allowed")
  }

  W
}
