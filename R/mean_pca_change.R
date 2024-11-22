
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
#' @details This functions performs structural break analysis for the functional
#'  data using an fPCA based initial dimension reduction. It is recommended
#'  that the dimension of the subspace, \code{d}, that the functional
#'  observations are projected onto should be selected based on TVE using
#'  \code{\link{pick_dim}}.
#'
#' @references Berkes, I., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Detecting changes in the mean of functional observations}
#'  Journal of the Royal Statistical Society, Series B 71, 927–946
#' @references Aue, A., Gabrys, R.,Horvath, L. & P. Kokoszka (2009).,
#'  \emph{Estimation of a change-point in the mean function of functional data}
#'  Journal of Multivariate Analysis 100, 2254–2269.
#'
#' @examples
#' res <- mean_pca_change(generate_brownian_bridge(200,seq(0,1,length.out=10)))
#' res1 <- mean_pca_change(generate_ou(20,200))
#' res2 <- mean_pca_change(funts(electricity))
mean_pca_change <- function(X, TVE=0.95,
                            M=1000, h=0, K=bartlett_kernel){
  X <- .check_data(X)
  X <- center(X)

  D <- nrow(X$data)
  n <- ncol(X$data)

  pca_X <- pca(X, TVE=TVE)

  d <- length(pca_X$sdev)

  eta.hat <- as.matrix(pca_X$x)

  ## Test Statistic
  Snd_tmp <- rep(NA,d)
  for(l in 1:d){
    inner_sums <- rep(0, n)
    for(k in 1:n){
      inner_sums[k] <- sum(eta.hat[1:k,l]) - k/n * sum(eta.hat[,l])
    }
    Snd_tmp[l] <- 1/pca_X$sdev[l]^2 * sum(inner_sums^2)
  }

  Snd <- 1/n^2 * sum(Snd_tmp)

  ## Pick k*
  # Sigma.hat = .long_run_cov(X, h=h,K = K)[1:d,1:d]
  kappa <-
    sapply(1:n,function(k, eta.hat, n){
      colSums(eta.hat[1:k, , drop=FALSE]) - k/n * colSums(eta.hat)
    },eta.hat=eta.hat, n=n)
  Q_nk <- rep(NA,n)
  for(k in 1:n){
    Q_nk[k] <- 1/n * ( t(kappa[,k]) %*% diag(1/pca_X$sdev^2) %*% kappa[,k] )
  }
  k_star <- which.max(Q_nk)

  # ASMPYD
  values <- sapply(1:M, function(k,d,n){
    B.Bridges <- rep(0, d)
    for(j in 1:d){
      B.Bridges[j] <- dot_integrate( sde::BBridge(0,0,0,1,n-1)^2 )
    }
    sum(B.Bridges)
  }, d=d,n=1000)

  pvalue <- sum(Snd <= values) / M

  dat.b <- funts(X$data[,1:k_star])
  dat.a <- funts(X$data[,(k_star+1):n])
  mean.b <- mean(dat.b)
  mean.a <- mean(dat.a)
  delta <- mean.a - mean.b

  ## Plots
  plot1 <- rainbow_plot(X$data,CPs=k_star)
  # distribution_plot(X,CPs=k_star) ## TODO:: color bands

  plot2 <-
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x=X$intraobs, y=delta)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype='dotted') +
    ggplot2::theme_bw() +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  list('data_plot'=plot1, 'mean_diff_plot'=plot2,
       pvalue = pvalue, change = k_star,
       DataBefore = dat.b, DataAfter = dat.a,
       change_fun = delta)
}



#' Compute S_n Test statistic for PCA change
#'
#' @param eta.hat TODO
#' @param k TODO
#' @param n TODO
#'
#' @return Numeric test statistic
#'
#' @keywords internal
#' @noRd
.S_n_pca <- function(eta.hat, k, n){
  # TODO:: Add
  # normalizer = ((k/n) * ((n-k) / n))^(-0.5)
  if(is.null(dim(eta.hat))){
    eta.bar <- sum(eta.hat)/n
    out <- sum(eta.hat[1:k]) - k*eta.bar
  } else {
    eta.bar <- as.matrix(rowSums(eta.hat)/n)
    out <- rowSums(as.matrix(eta.hat[, 1:k])) - k*eta.bar
  }

  out #* normalizer
}

#' Asymptotic Threshold for mean PCA change
#'
#' @param N description
#' @param d description
#'
#' @return Numeric asymptotic threshold
#'
#' @keywords internal
#' @noRd
.asymp_pca <- function(N, d){
  B.Bridges <- matrix(0,nrow = d,ncol = N)
  for(j in (1:d)){
    B.Bridges[j,] <- sde::BBridge(0,0,0,1,N-1)^2
  }
  max(colSums(B.Bridges))
}
