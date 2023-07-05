#' Scalar Characteristic Function Change Point Analysis
#'
#' @param Y
#' @param gam
#' @param nSims
#' @param alpha_val
#'
#' @return
#' @export
#'
#' @examples
scalarDetection <- function(Y, gam = 0.5, nSims = 200,
                             alpha_val=0.05){
  n <- length(Y)

  Tn_s <- rep(NA,n-1)
  for(k in 1:(n-1)){
    Tn_s[k] <- ((k*(n-k))/n^2)^gam * ((k*(n-k))/n) *
      integrate(function(t,Y,k){abs(.phi_k(Y,t,k)-.phi_k0(Y,t,k))^2 * .w(t)},
                lower=0, upper=1, Y=Y, k=k)[[1]]
  }
  Tn <- max(Tn_s)

  Tn_permute <- sapply(rep(NA,nSims),
                       function(tmp, n, gam, Y){

                         Tn_tmp_s <- sapply(1:(n-1), function(k,n,gam,Y){
                           ((k*(n-k))/n^2)^gam * ((k*(n-k))/n) *
                             integrate(function(t,Y,k){abs(.phi_k(Y,t,k)-.phi_k0(Y,t,k))^2 * .w(t)},
                                       lower=0, upper=1, Y=Y[sample(1:length(Y))],k=k)[[1]]
                         },n=n,gam=gam,Y=Y)

                         max(Tn_tmp_s)
                       }, n=n, gam=gam, Y=Y)

  list('dat'=Y,
       'Tn'=Tn,
       'Tn_permute'=Tn_permute,
       'pval'=1 - ecdf(Tn_permute)(Tn),
       'Loc'=which.max(Tn_s),
       'loc_if_pval'=ifelse(1 - ecdf(Tn_permute)(Tn) < alpha_val,
                            which.max(Tn_s), NA))
}


.phi_k <- function(Y, t, k){
  sumVal <- 0
  for(j in 1:k){
    sumVal <- sumVal +
      exp(complex(real=0,imaginary = 1) * t * Y[j])
  }

  1/k * sumVal
}

.phi_k0 <- function(Y, t, k){
  n <- length(Y)
  sumVal <- 0

  for(j in (k+1):n){
    sumVal <- sumVal +
      exp(complex(real=0,imaginary = 1) * t * Y[j])
  }

  1/(n-k) * sumVal
}

.w <- function(t, a=1){exp(-a*abs(t))}
