#' Change Point Confidence Intervals
#'
#' @param X
#' @param CPs
#' @param K Default is bartlett_kernel.
#' @param h Default is 3.
#' @param kappa Default is 0.5.
#' @param M Default is 1000.
#' @param alpha Default is 0.05
#' @param method Options are 'distribution' and 'simulation'. Default is distribution.
#'
#' @references Horváth, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (First edition.). Springer. \url{https://doi.org/10.1007/978-3-031-51609-2}
#'
#' @references Aue, A., Rice, G., & Sönmez, O. (2018). Detecting and dating structural
#'  breaks in functional data without dimension reduction. Journal of the Royal
#'  Statistical Society. Series B, Statistical Methodology, 80(3), 509–529.
#'  \url{https://doi.org/10.1111/rssb.12257}
#'
#' @return
#' @export
#'
#' @examples
#' X <- cbind(generate_brownian_motion(100,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+1000)
#' compute_confidence_interval(X,CPs = 100)
#' compute_confidence_interval(X,CPs=100,method = 'simulation')
#'
#' X <- cbind(generate_brownian_motion(100,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.5)
#' compute_confidence_interval(X,100,alpha = 0.1)
#' compute_confidence_interval(X,CPs=100,alpha = 0.1,method = 'simulation')
#'
#' X <- generate_brownian_motion(200,v=seq(0,1,0.05))
#' compute_confidence_interval(X,100)
#' compute_confidence_interval(X,100,method = 'simulation')
#'
#' X <- cbind(generate_brownian_motion(200,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.1,
#'            generate_brownian_motion(150,v=seq(0,1,0.05))$data-0.05)
#' compute_confidence_interval(X,c(200,300))
#'
#' # set.seed(12345)
#' # bs <- binary_segmentation(electricity,statistic = 'Tn',method = 'Boot')
#' compute_confidence_interval(X = electricity, CPs = c(66, 144, 204, 243, 305),alpha = 0.1)
compute_confidence_interval <- function(X, CPs, K=bartlett_kernel,
                                h=3, kappa=0.5, M=5000,
                                alpha=0.1, method='distribution'){
  method <- .verify_input(method, c('distribution','simulation'))

  X <- .check_data(X)

  CPs_ext <- c(0,CPs,ncol(X$data))
  r <- nrow(X$data)
  results <- data.frame('change' = rep(NA,length(CPs)),
                        'lower' =  NA, 'upper' = NA)

  for(i in 1:length(CPs)){
    CP <- CPs[i]
    CP_simple <- CP-CPs_ext[i]
    idx1 <- (CPs_ext[i]+1):CP
    idx2 <- (CP+1):CPs_ext[i+2]
    N <- length( c(idx1,idx2) )

    X_k1 <- rowMeans(X$data[,idx1])
    X_k2 <- rowMeans(X$data[,idx2])

    e <- X_k2 - X_k1
    e_l2norm <- dot_integrate_uneven(e^2, r = X$intraobs)
    eps_hat <- matrix(ncol=N, nrow=r)
    eps_hat[,1:CP_simple] <- X$data[,idx1] - X_k1
    eps_hat[,(CP_simple+1):N] <- X$data[,idx2] - X_k2

    g <- rep(NA,N)
    for(j in 1:N){
      g[j] <- dot_integrate_uneven(
        eps_hat[,j] * e / sqrt(e_l2norm),
        r = X$intraobs)
    }

    ## LRV
    iters <- (1-N):(N-1)
    Kvals <- K(iters/h)
    data_tmp <- data.frame('K'=Kvals[Kvals>0],
                           'l'=iters[which(Kvals>0)])
    values <- apply(data_tmp, MARGIN = 1,
                    function(ker_loc, g, N) {
                      if (ker_loc[2] >= 0) {
                        rs <- 1:(N - ker_loc[2])
                        con_gamma_l <- g[rs]*g[rs+ker_loc[2]]
                      } else {
                        rs <- (1 - ker_loc[2]):N
                        con_gamma_l <- g[rs]*g[rs+ker_loc[2]]
                      }

                      sum(ker_loc[1] * con_gamma_l) * 1 / (N-abs(ker_loc[2]))
                    },
                    g = g, N = N
    )
    tau2 <- sum(values)

    brown_v <- seq(0,20,0.1)
    double_t <- c(rev(-brown_v[-1]), brown_v)
    theta <- CP_simple / N

    if(method=='simulation'){

      m <- rep(0, length(double_t))
      for(t_idx in 1:length(m)){
        if(double_t[t_idx] < 0){
          m[t_idx] <- (1 - kappa)*(1 - theta) + kappa*theta
        }else if(double_t[t_idx] > 0){
          m[t_idx] <- (1 - kappa)*theta + kappa*(1 - theta)
        }
      }

      # xi <- rep(NA,M)
      # for(j in 1:M){
      #   W <- c(rev(generate_brownian_motion(1, v=brown_v)$data[-1,1]),
      #          generate_brownian_motion(1, v=brown_v)$data[,1])
      #
      #   xi[j] <- double_t[which.max(W - abs(double_t) * m)]
      # }

      W1 <- generate_brownian_motion(M, v=brown_v)$data[length(brown_v):2,]
      W2 <- generate_brownian_motion(M, v=brown_v)$data
      W <- rbind(W1,W2)
      xi <- apply(W,MARGIN = 2,function(w,double_t,m){
        double_t[which.max(w - abs(double_t) * m)]
      },double_t=double_t,m=m)

      thetas0 <- quantile(xi, probs = c(alpha/2, 1-alpha/2))

    } else if(method=='distribution'){
      # Pg 40 (pdf 52), information
      .h_func <- function(u,x,y){
        2*x * (x+2*y) * (1- pnorm( (x+2*y)*sqrt(u) )) *
          exp(2*y * (x+y) * u) - 2*x^2 * (1-pnorm( x*sqrt(u) ))
      }
      .f_distribution <- function(t, kappa, theta){
        ifelse(t<=0,
          suppressWarnings(.h_func(-t, (1-kappa)*(1-theta)+kappa*theta, (1-kappa)*theta + kappa*(1-theta))),
          suppressWarnings(.h_func(t, (1-kappa)*theta+kappa*(1-theta), (1-kappa)*(1-theta) + kappa*theta))
        )
      }

      thetas0 <- rep(NA,2)
      f_CDF <- function(up, kappa, theta, alpha){
        abs(integrate(.f_distribution,-20,up, kappa=kappa, theta=theta)$value - alpha)
      }
      # thetas0[1] <- -10.63
      #   # dot_integrate_uneven(.f_distribution(seq(-25,-10.63,0.01),kappa,theta),
      #   #                      seq(-25,-10.63,0.01))
      # thetas0[2] <- 11.4
      #   # dot_integrate_uneven(.f_distribution(seq(-25,11.4,0.01),kappa,theta),
      #   #                      seq(-25,11.4,0.01))

      tmp <- optimize(f_CDF, interval = c(-20,20),
                      kappa=kappa, theta=theta, alpha=alpha/2)$minimum
      thetas0[1] <- tmp#.f_distribution(tmp,kappa,theta)
      tmp <- optimize(f_CDF, interval = c(-20,20),
                      kappa=kappa, theta=theta, alpha=1-alpha/2)$minimum
      thetas0[2] <- tmp#.f_distribution(tmp,kappa,theta)

    }else{
      stop('Method incorrectly specified',call. = FALSE)
    }

    results[i,] <- c(CP,
                     as.numeric( CP + tau2*thetas0[1] / e_l2norm ),
                     as.numeric( CP + tau2*thetas0[2] / e_l2norm ) )

  }

  results
}

