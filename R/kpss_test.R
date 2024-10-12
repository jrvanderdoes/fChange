#' Functional KPSS Test
#'
#' @param X description
#' @param method description
#' @param boot_method description
#' @param M description
#' @param h description
#' @param TVE description
#' @param replace description
#'
#' @return
#' @export
#'
#' @references Chen, Y., & Pun, C. S. (2019). A bootstrap-based KPSS test for
#'  functional time series. Journal of Multivariate Analysis, 174, 104535.
#' @references Kokoszka, P., & Young, G. (2016). KPSS test for functional time
#'  series. Statistics, 50(5), 957-973.
#'
#' @examples
#' compute_kpss(generate_brownian_motion(100,v=seq(0,1,length.out=20)))
#' compute_kpss(generate_brownian_motion(100,v=seq(0,1,length.out=20)),
#'              method="BS")
#' compute_kpss(generate_brownian_motion(100,v=seq(0,1,length.out=20)),
#'              method="BS", boot_method='overlapping')
compute_kpss <-function(X, method="MC", boot_method='seperate',
                        M=1000, h = 3, TVE=0.95, replace=TRUE){
  X <- .check_data(X)

  N <- ncol(X$data)
  r <- nrow(X$data)

  tmp <- .compute_fpss_test_stat(X)
  RN <- tmp$test_statistic
  hat_eta <- tmp$hat_eta

  if (method=="MC"){

    pca_X <- pca(X, TVE=TVE)

    sim_RNs <- sapply(1:M,function(m,eigs,v){
      sum(eigs *
            dot_integrate_col(.generate_second_BB(length(eigs), v = v)$data^2) )
    },eigs=pca_X$sdev^2, v=X$intraobs)

  } else if(method=='BS'){
    xi_hat <- (12/(N*(N^2-1))) *
      rowSums(matrix(rep(1:N-(N+1)/2, r ),
             nrow=r, ncol=N, byrow = TRUE ) * X$data)
    mu_hat <- rowMeans(X$data) - xi_hat * ((N+1)/2)

    boot_etas <- .bootstrap(hat_eta, blockSize=h, M=M,
                           type=boot_method, replace=replace)
    boot_X <- list()
    sim_RNs <- rep(NA, M)
    for(i in 1:M){
      boot_X[[i]] <- matrix(mu_hat,nrow = r,ncol = N) +
        matrix(1:N,ncol=N, nrow=r,byrow = TRUE) * xi_hat + boot_etas[[i]]

      tmp <- .compute_fpss_test_stat(funts(boot_X[[i]]))
      sim_RNs[i] <- tmp$test_statistic
    }

  }

  list(statistic = RN, pvalue = sum(RN <= sim_RNs)/M)
}



#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
#'
#' @keywords internal
.compute_fpss_test_stat <- function(X){
  N <- ncol(X$data)
  r <- nrow(X$data)

  iX <- rowSums(t( 1:N * t(X$data) ))

  hat_eta <- X$data +
    t( ( (1:N-(N+1)/2) * 6/(N-1) - 1 )  %*% t(mean(X))) -
    t( ( 12/(N*(N^2-1))*(1:N-(N+1)/2) ) %*% t(iX) )

  ZN <- 1/sqrt(N) * cumsum(funts(hat_eta))$data
  # TODO: UNeven integregrate
  RN <- dot_integrate(dot_integrate_col(ZN^2)) * r

  list('test_statistic'=RN,'hat_eta'=hat_eta)
}

#' Title
#'
#' @param M
#' @param v
#'
#' @return
#' @export
#'
#' @examples
.generate_second_BB <- function(M,v){
  # TODO:: Add to generation script
  W <- generate_brownian_motion(M,v=v)
  V <- W$data +
    (2*v-3*v^2) %*% t(W$data[nrow(W$data),]) +
    (-6*v + 6*v^2) %*% t(dot_integrate_col(W$data))

  funts(X=V,intraobs = v)
}
