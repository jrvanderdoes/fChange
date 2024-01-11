generate_data_n <- function(N=100, resolution=30){

  X <- data.frame(matrix(nrow=resolution,ncol=N))
  W <- computeSpaceMeasuringVectors(N,space = 'BM',X=X)

  K <- function(t,s,c){c*exp(-(t^2+s^2)/2)}

  dot_integrate_col()
  for(i in 1:N){
    X[,i] <- dot_integrate(K(,s)* X[,i-1] + W[,i])
  }
}
