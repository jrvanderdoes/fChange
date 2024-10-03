# Convert X to pca
if(is.null(CPs)){
  X_pca <- pca(X,TVE = TVE)
  D <- min(ncol(X_pca$x), max.d)
}else{
  CPs <- unique(c(0, CPs, ncol(X)))
  # max_n <- max((CPs-lag(CPs))[-1])
  # D <- length(CPs)-1
  # dat <- data.frame(matrix(nrow=max_n,ncol=D))
  X_demean <- X
  for(d in 1:(length(CPs)-1)){
    X_tmp <- funts(X$data[,(CPs[d]+1):CPs[d+1]],intraobs = X$intraobs)
    X_demean$data[,(CPs[d]+1):CPs[d+1]] <- X_tmp$data - mean(X_tmp)
  }
  X_pca <- pca(X_demean, TVE = TVE)
  D <- min(ncol(X_pca$x), max.d)
}

dat <- X_pca$x[,1:D]
