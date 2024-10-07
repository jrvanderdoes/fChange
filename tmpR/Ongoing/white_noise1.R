.compute_sqrtMat <- function(X,W,J,h,K){
  x <- seq(0,1,length.out=J)

  # Compute Gamma Matrix
  covMat <- .estimCovMat(X, W, x, h, K)

  # TODO:: Test "GauPro::sqrt_matrix()" which is slightly faster
  tryCatch({
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  }, error = function(e){
    covMat <- round(covMat, 15)
    covMat_svd <- La.svd(covMat)
    sqrtD <- diag(sqrt(covMat_svd$d))
  })

  Rfast::mat.mult(
    Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
  )
}

.estimCovMat <- function(X, W, x,
                         h = 3, K = bartlett_kernel) {
  M <- ncol(W)

  D11 <- D12_tD21 <- D22 <- data.frame(matrix(ncol = M, nrow = M))
  funs <- c(cos, sin)
  for (i in 1:M) {
    for (j in i:M) {
      D11[i, j] <- D11[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[1][[1]], vp = W[, j]
        )
      D12_tD21[i, j] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
        )
      if (j > i) {
        D12_tD21[j, i] <- .estimD(
          K = K, h = h, X = X,
          lfun = funs[1][[1]], v = W[, j],
          lfunp = funs[2][[1]], vp = W[, i]
        )
      }

      D22[i, j] <- D22[j, i] <-
        .estimD(
          K = K, h = h, X = X,
          lfun = funs[2][[1]], v = W[, i],
          lfunp = funs[2][[1]], vp = W[, j]
        )
    }
  }

  tmp1 <- x %*% t(x)
  minDF <- matrix(1, nrow = length(x), ncol = length(x))
  for (i in 1:length(x)) {
    for (j in 1:length(x)) {
      minDF[i, j] <- (min(x[i], x[j]) - tmp1[i, j])
    }
  }

  fastmatrix::kronecker.prod(
    as.matrix(minDF),
    rbind(
      cbind(as.matrix(D11), as.matrix(D12_tD21)),
      cbind(as.matrix(t(D12_tD21)), as.matrix(D22))
    )
  )
}
