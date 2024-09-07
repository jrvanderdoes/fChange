
#' Change Point Analysis - Dimension Reduction
#'
#' This function tests whether there is a significant change in the mean
#'  function of functional data, and it gives an estimate of the location
#'  of the change. The procedure will reduce the dimension of the
#'  functional data using functional principal component analysis and will
#'  use the leading principal curves which explain \code{TVE} total variance
#'  to carry out the change point analysis.
#'
#' @param fdobj A functional data object of class '\code{fd}'
#' @param d Number of principle components
#' @param M Number of monte carlo simulations to get the critical values. The default value is \code{M=1000}
#' @param  h The window parameter for the estimation of the long run covariance kernel. The default
#' value is \code{h=0}, i.e., it assumes iid data
#' @param plot If \code{TRUE} plot of the functional data before and after the estimated change and plot of the
#' estimated change function is given
#' @param ... Further arguments to pass
#'
#' @export
#'
#' @return
#'\item{\code{pvalue}}{
#' An approximate p value for testing whetehr there is a significant change in the mean function
#'}
#'\item{\code{change}}{
#' Estimated change location
#'}
#'\item{\code{DataBefore}}{
#' Data before the estimated change
#'}
#'\item{\code{DataAfter}}{
#' Data after the estimated change
#'}
#'\item{\code{MeanBefore}}{
#' Mean function before the estimated change
#'}
#'\item{\code{MeanAfter}}{
#' Mean function after the estimated change
#'}
#'\item{\code{change_fun}}{
#' Estimated change function
#'}
#'
#' @seealso \code{\link{change_FF}}
#'
#'@details This functions performs structural break analysis for the functional data using an fPCA based initial dimension reduction. It
#'is recommended that the dimension of the subspace, \code{d}, that the functional observations are projected onto should be selected based on
#'TVE using \code{\link{pick_dim}}.
#'
#'@references Berkes, I., Gabrys, R.,Hovarth, L. & P. Kokoszka (2009)., \emph{Detecting changes in the mean of functional observations}
#' Journal of the Royal Statistical Society, Series B 71, 927–946
#' @references Aue, A., Gabrys, R.,Hovarth, L. & P. Kokoszka (2009)., \emph{Estimation of a change-point in the mean function
#' of functional data}Journal of Multivariate Analysis 100, 2254–2269.
#'
#' @examples
#' mean_pca_change(funts(electricity))
mean_pca_change <- function(X, TVE=0.95, M=1000, h=0, K=bartlett_kernel){
  stop('Errors still remaining')
  fdata <- fun_IID(n=200, nbasis=21)
  X <- fda::eval.fd(seq(0,1,length.out=21),fdata)
  X <- electricity
  ## TODO:: FIX EXAMPLE
  X <- .check_data(X)
  X <- center(X)

  # D <- nrow(X$data)
  n <- ncol(X$data)

  pca_X <- pca(X, TVE)
  d <- length(pca_X$sdev)

  eta.hat <- t(pca_X$rotation)

  Sigma.hat <- .long_run_var(X$data, h=h, K=K)[1:d, 1:d]
  TT <- sapply(1:n,
               function(k,n,eta.hat) {
                 1/n *( t(.S_n_pca(eta.hat = eta.hat,k = k,n = n)) %*%
                          solve(Sigma.hat) %*%
                          .S_n_pca(eta.hat = eta.hat,k = k,n = n))},
               n=n, eta.hat=eta.hat)
  Tn <- max(TT)
  k.star <- min(which(TT==max(TT)))

  Values <- sapply(1:M, function(k,d) {.asymp_pca(n, d)},d=d)
  # z = Tn<=Values
  # p = length(z[z==TRUE])/length(z)
  p <- sum(Tn <= Values) / M # Compute p-value
  p
  Tn
  Values
  # p <- 1-ecdf(Values)(Tn)
  dat.b <- funts(X$data[,1:k.star])
  dat.a <- funts(X$data[,(k.star+1):n])
  mean.b <- mean(dat.b)
  mean.a <- mean(dat.a)
  delta <- mean.a - mean.b

  # par(mfrow=c(1,2))
  ## add before/after
  plot1 <- .plot_stack(X,CPs=k.star) ## TODO:: lower, add mean
  # .plot_substack(X,CPs=k.star) ## TODO:: color bands

  plot2 <- plot(delta, main="Estimated Change Function",
                ylab="values", type='l')
  list('data_plot'=plot1, 'mean_diff_plot'=plot2,
       pvalue = p , change = k.star,
       DataBefore = dat.b, DataAfter = dat.a,
       change_fun = delta)
}



#' Title
#'
#' @param k
#'
#' @return
.S_n_pca = function(eta.hat, k, n){
  # TODO:: Add
  # normalizer = ((k/n) * ((n-k) / n))^(-0.5)
  if(d==1){
    eta.bar = sum(eta.hat)/n
    out = sum(eta.hat[1:k]) - k*eta.bar
  } else {
    eta.bar = as.matrix(rowSums(eta.hat)/n)
    out = rowSums(as.matrix(eta.hat[, 1:k])) - k*eta.bar
  }

  out #* normalizer
}

#' Title
#'
#' @param N
#'
#' @return
.asymp_pca <- function(N, d){
  B.Bridges= matrix(0,d,(N))
  for(j in (1:d)){
    B.Bridges[j,] <- sde::BBridge(0,0,0,1,N-1)^2
  }
  max(colSums(B.Bridges))
}
