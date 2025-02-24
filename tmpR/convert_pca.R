# Convert X to pca
if(is.null(changes)){
  X_pca <- pca(X,TVE = TVE)
  D <- min(ncol(X_pca$x), max.d)
}else{
  changes <- unique(c(0, changes, ncol(X)))
  # max_n <- max((changes-lag(changes))[-1])
  # D <- length(changes)-1
  # dat <- data.frame(matrix(nrow=max_n,ncol=D))
  X_demean <- X
  for(d in 1:(length(changes)-1)){
    X_tmp <- funts(X$data[,(changes[d]+1):changes[d+1]],intratime = X$intratime)
    X_demean$data[,(changes[d]+1):changes[d+1]] <- X_tmp$data - mean(X_tmp)
  }
  X_pca <- pca(X_demean, TVE = TVE)
  D <- min(ncol(X_pca$x), max.d)
}

dat <- X_pca$x[,1:D]
