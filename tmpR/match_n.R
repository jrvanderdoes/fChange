
detect_changepoint_singleCov_n <- function(X, nSims = 2000, x = seq(0, 1, length.out = 30),
                                           h = 3,
                                           K = bartlett_kernel,
                                           space = "BM",
                                           silent = FALSE, TN_M = 100000, Cov_M = 40,
                                           sqrtMat=NULL) {
  # Determine Number of Iterations
  val_Tn <- compute_Tn_n(X, length(x), M = TN_M, space = space)
  MJ <- Cov_M * length(x)

  if(is.null(sqrtMat)){
    sqrtMat <- .compute_sqrtMat(Cov_M=Cov_M, space=space, X=X,x=x,h=h,K=K)
    # # Generate Noise
    # W <- computeSpaceMeasuringVectors(Cov_M, space, X)
    #
    # # Variables
    # MJ <- Cov_M * length(x)
    #
    # # Compute Gamma Matrix
    # # Make sure x is resolution of data
    # # TODO:: Try resolution for n or clump tn to map same
    # covMat <- .estimCovMat(X, x, Cov_M, h, K, W)
    # covMat <- round(covMat, 15)
    # covMat_svd <- La.svd(covMat)
    # sqrtD <- diag(sqrt(covMat_svd$d))
    # sqrtD[is.na(sqrtD)] <- 0
    # #sqrtMat <- eigenMapMatMult(
    # #  eigenMapMatMult(covMat_svd$u, sqrtD),t(covMat_svd$v))
    # sqrtMat <- Rfast::mat.mult(
    #   Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
    # )
    # rm(W, covMat, covMat_svd, sqrtD)
  }

  gamProcess <- c()
  nIters <- nSims / 100
  gamProcess <- sapply(1:nIters, FUN = function(tmp, MJ, sqrtMat) {
    # (After trans + mult) Rows are iid MNV
    # mvnorms <- sapply(1:100, function(m, x) {
    #   stats::rnorm(x)
    # }, x = 2 * MJ)
    mvnorms <- matrix(stats::rnorm(100*2*MJ),ncol=100,nrow=2*MJ)
    mvnorms <- Rfast::mat.mult(sqrtMat, mvnorms)

    gamVals <- mvnorms[1:MJ, ] + complex(imaginary = 1) * mvnorms[MJ + 1:MJ, ]
    # gamVals <- mvnorms[(M+1):MJ,] + complex(imaginary = 1)*mvnorms[MJ+1:(MJ-M),]

    # Estimate value
    apply(abs(t(gamVals))^2, MARGIN = 1, dot_integrate)
  }, MJ = MJ, sqrtMat = sqrtMat)
  gamProcess <- as.vector(gamProcess)

  list(
    "pval" = 1 - stats::ecdf(gamProcess)(val_Tn),
    "gamProcess" = gamProcess,
    "value" = val_Tn
  )
}



#############################################



compute_Tn_n <- function(X, n, M = 100000, W = NULL, space = "BM", ...) {
  n <- ncol(X)

  if (is.null(W)) {
    W <- computeSpaceMeasuringVectors(M = M, X = X, space = space)
  } else {
    M <- ncol(W)
  }


  Zn <- .Zn_n(W, X, .select_n(1:ncol(X),n))
  intVal <- dot_integrate_col(abs(Zn)^2)

  sum(intVal)/M
}


.Zn_n <- function(W, X, ns) {
  fhat_vals <- as.matrix(.fhat_all(X, W))

  sqrt(length(ns)) * (fhat_vals[ns,] - (ns/length(ns)) %o% fhat_vals[nrow(fhat_vals),])
}

.select_n <- function(vals,n){
  vals[round(seq(1,length(vals),length.out=n))]
}
