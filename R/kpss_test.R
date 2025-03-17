#' Functional KPSS Test
#'
#' Compute the KPSS statistic for functional data.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param method String (default MC) for the method: Monte Carlo simulation
#'  (\code{MC}), block resampling (\code{block}), and sliding window
#'  permutation (\code{sliding})
#' @param M Number of simulations (default 1000) to estimate theoretical distribution
#' @param h Numeric (default 3) for the block size when using a resample test
#' @param TVE Numeric (default 1) for selecting the number of principle
#'  components / eigenvalues
#' @param replace Boolean (default TRUE) to indicate if blocks should be
#'  selected with replacement when using a resample test
#'
#' @return List with statistic, p-value and theoretical simulations, respectively
#' @export
#'
#' @references Chen, Y., & Pun, C. S. (2019). A bootstrap-based KPSS test for
#'  functional time series. Journal of Multivariate Analysis, 174, 104535.
#' @references Kokoszka, P., & Young, G. (2016). KPSS test for functional time
#'  series. Statistics, 50(5), 957-973.
#'
#' @examples
#' kpss_test(generate_brownian_motion(100,v=seq(0,1,length.out=20)))
#' kpss_test(generate_brownian_motion(100,v=seq(0,1,length.out=20)),
#'              method="block")
#' kpss_test(generate_brownian_motion(100,v=seq(0,1,length.out=20)),
#'              method='sliding')
kpss_test <-function(X, method=c("MC",'block','sliding'),
                        M=1000, h = 3, TVE=1, replace=TRUE){
  X <- dfts(X)
  method <- .verify_input(method, c("MC",'block','sliding'))

  N <- ncol(X$data)
  r <- nrow(X$data)

  tmp <- .compute_kpss_test_stat(X)
  RN <- tmp$test_statistic
  hat_eta <- tmp$hat_eta

  if (method=="MC"){

    # pca_X <- pca(X, TVE=1)$sdev^2
    # eigs <- est_eigenvalue(t(hat_eta), K, 2/5)
    eigs <- pca(dfts(hat_eta), TVE=TVE)$sdev^2

    sim_RNs <- sapply(1:M,function(m,eigs,v){
      sum(eigs *
            dot_integrate_col(
              .generate_second_level_brownian_bridge(length(eigs), v = v)$data^2) )
    },eigs=eigs, v=X$fparam)

  } else if(method=='block' || method=='sliding'){
    xi_hat <- (12/(N*(N^2-1))) *
      rowSums(matrix(rep(1:N-(N+1)/2, r ),
             nrow=r, ncol=N, byrow = TRUE ) * X$data)
    mu_hat <- rowMeans(X$data) - xi_hat * ((N+1)/2)

    boot_etas <- .bootstrap(hat_eta, blocksize=h, M=M,
                           type=ifelse(method=='block','separate','overlapping'),
                           replace=replace)
    boot_X <- list()
    sim_RNs <- rep(NA, M)
    for(i in 1:M){
      boot_X[[i]] <- matrix(mu_hat,nrow = r,ncol = N) +
        matrix(1:N,ncol=N, nrow=r,byrow = TRUE) * xi_hat + boot_etas[[i]]

      tmp <- .compute_kpss_test_stat(dfts(boot_X[[i]]))
      sim_RNs[i] <- tmp$test_statistic
    }

  }

  list(statistic = RN, pvalue = sum(RN <= sim_RNs)/M, sims=sim_RNs)
}



#' Compute KPSS Test statistic
#'
#' @inheritParams kpss_test
#'
#' @return List with test statistics and non-integrated values
#'
#' @keywords internal
#' @noRd
.compute_kpss_test_stat <- function(X){
  N <- ncol(X$data)
  r <- nrow(X$data)

  iX <- rowSums(t( 1:N * t(X$data) ))

  hat_eta <- X$data +
    t( ( (1:N-(N+1)/2) * 6/(N-1) - 1 )  %*% t(mean(X))) -
    t( ( 12/(N*(N^2-1))*(1:N-(N+1)/2) ) %*% t(iX) )

  # ZN <- 1/sqrt(N) * cumsum(dfts(hat_eta))$data
  ZN <- .kpss_partial_sum(hat_eta)
  RN <- dot_integrate(dot_integrate_col(ZN^2, X$fparam))# * r

  list('test_statistic'=RN,'hat_eta'=hat_eta)
}


#' KPSS Partial Sum Function
#'
#' Function for partial sum in KPSS test
#'
#' @inheritParams kpss_test
#'
#' @return Data.frame of eta_hat
#'
#' @keywords internal
#' @noRd
.kpss_partial_sum <- function(X){
  tX <- t(X)
  N <- nrow(tX)
  res <- ncol(tX)
  teta <- matrix(rep(0,res*res),ncol=res)

  if (floor(N/res)>0) {
    for (i in 1:floor(N/res)){
      teta[1,] <- teta[1,] + tX[i,]/sqrt(N)
    }
  }

  if (res>1){
    for (i in 2:res){
      teta[i,] <- teta[i-1,]
      if (floor(N*(i-1)/res) < floor(N*i/res)) {
        for (j in (floor(N*(i-1)/res)+1):floor(N*i/res)){
          teta[i,] <- teta[i,] + tX[j,]/sqrt(N)
        }
      }
    }
  }

  t(teta)
}


# est_eigenvalue<-function(este,K,h_power=2/5){
#   N <- nrow(este)
#   res <- ncol(este)
#
#   gamma<-array(0, dim=c(N,res,res))
#   for (i in 1:N){
#     for (j in i:N){
#       gamma[i,,] <- gamma[i,,] + ( este[j,] %*% t(este[j-i+1,]) )/N
#     }
#   }
#
#   cc <- gamma[1,,]
#   for (i in 2:N){
#     cc <- cc + K( (i-1)/(N^(h_power)) ) * (gamma[i,,]+t(gamma[i,,]))
#   }
#
#   eigen(cc/res)$values
# }
