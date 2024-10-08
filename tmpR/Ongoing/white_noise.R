white_noise_test <- function(X, M = 20, J=50, l=1, nSims = 1000, h = 3,
                             K = bartlett_kernel, space = "BM", silent = FALSE) {
  # Generate Noise
  W1 <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  W2 <- computeSpaceMeasuringVectors(M = M, X = X, space = space)

  # Determine Number of Iterations
  double_integrate <- .wn_test_stat(X = X, W1 = W1, W2 = W2, l=l, J = J)

  MJ <- M * J

  X_W1W2 <- .compute_sqrtMat(X, W1, W2, J, h, K)
  X_W1 <- .compute_sqrtMat(X$data[,-ncol(X)],W1, J, h, K)
  X_W2 <- .compute_sqrtMat(X$data[,-1],W2, J, h, K)

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat, M, J) {
    # (After trans + mult) Rows are iid MNV
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms1 <- Rfast::mat.mult(t(mvnorms), sqrtMat)
    # 0 0 0 0 (M times) ... .... 1 1 1 1 (M times)

    gamVals <- mvnorms1[,1:MJ] + complex(imaginary = 1) * mvnorms1[,MJ + 1:MJ]

    # Integrate out M
    results <- matrix(NA,nrow=100, ncol=J)
    for(i in 1:J){
      results[,i] <- apply(abs(gamVals[,(i-1)*M +1:M])^2,
                           MARGIN=1, mean)
    }

    # Estimate value
    list(apply(results, MARGIN = 1, dot_integrate),
         apply(results, MARGIN = 1, max) )
  }, MJ = MJ, sqrtMat = sqrtMat, M=M, J=J)
  gamProcessTn <- as.vector(unlist(gamProcess[1,]))
  gamProcessMn <- as.vector(unlist(gamProcess[2,]))

  list(
    'Tn'=list(
      "pval" = 1 - stats::ecdf(gamProcessTn)(val_Tn),
      "gamProcess" = gamProcessTn,
      "value" = val_Tn
    ),
    'Mn'=list(
      "pval" = 1 - stats::ecdf(gamProcessMn)(val_Mn$value),
      "gamProcess" = gamProcessMn,
      "value" = val_Mn$value
    ),
    'Location'=val_Mn$location
  )
}

.wn_test_stat <- function(X, W1, W2, l=1, J=50) {
  X <- funts(X)

  n <- ncol(X$data)

  Zn <- .Zn_final(W, X$data)
  ns <- .select_n(1:n, J)
  Zn <- Zn[ns,,drop=FALSE]
  # Integrate out W
  intVal <- rowMeans(abs(Zn)^2)#dot_integrate_col(t(abs(Zn)^2))

  # Integrate observations
  dot_integrate(intVal)
}

.Zn_final <- function(W, X) {
  n <- ncol(X)
  fhat_vals <- as.matrix(.fhat_all(X, W))

  sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])
}
