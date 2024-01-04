set.seed(1234)
change <- 0
data <- list()
nSims <- 1000
res <- 20
M <- 50

hVal <- 28

.generateFunctionalIID <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = FALSE)

  t(fiid)
}

data_tmp <- .generateFunctionalIID(res,15000)

# Generate Noise
W <- computeSpaceMeasuringVectors(M, 'BM', data_tmp)

# Variables
x <- seq(0, 1, length.out = res)
MJ <- M * length(x)

sqrtMat <- .compute_sqrtMat(Cov_M=M, space='BM', X=data_tmp,x=x,h=hVal,K=bartlett_kernel)
# # Compute Gamma Matrix
# covMat <- .estimCovMat(data_tmp, x, M, hVal, bartlett_kernel, W)
# covMat <- round(covMat, 15)
# covMat_svd <- La.svd(covMat)
# sqrtD <- diag(sqrt(covMat_svd$d))
# sqrtD[is.na(sqrtD)] <- 0
# sqrtMat <- Rfast::mat.mult(
#   Rfast::mat.mult(covMat_svd$u, sqrtD), covMat_svd$vt
# )
# rm(W, covMat, covMat_svd, sqrtD)

result <- data.frame('norm'=rep(NA,nSims),
                     'big'=NA,
                     'norm_nmatch'=NA,
                     'big_nmatch'=NA)

for(i in 1:nSims){
  cat(paste0(i,', '))
  set.seed(123+i*123)

  data[[i]] <- .generateFunctionalIID(res,15000)

  Cov_iid <- detect_changepoint_singleCov(
    X=data[[i]],x = x,TN_M = M,K = bartlett_kernel,
    nSims = 1000, h=hVal, space='BM', Cov_M=M, silent = T)
  Cov_iid1 <- detect_changepoint_singleCov(
    X=data[[i]], x = x,TN_M = M,K = bartlett_kernel,
    nSims = 1000, h=hVal, space='BM', Cov_M=M,
    silent = T, sqrtMat = sqrtMat)

  Cov_iidn <- detect_changepoint_singleCov_n(
    X=data[[i]], x = x,TN_M = M, K = bartlett_kernel,
    nSims = 1000, h=hVal, space='BM', Cov_M=M, silent = T)
  Cov_iidn1 <- detect_changepoint_singleCov_n(
    X=data[[i]], x = x,TN_M = M, K = bartlett_kernel,
    nSims = 1000, h=hVal, space='BM', Cov_M=M,
    silent = T,sqrtMat = sqrtMat)

  result[i,] <- c(Cov_iid$pval, Cov_iid1$pval,
                  Cov_iidn$pval, Cov_iidn1$pval)
  saveRDS(result,'C:/Users/jerem/Downloads/CE_check_new4.rds')
}

colMeans(result<=0.05)
