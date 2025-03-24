#' Change Point Confidence Intervals
#'
#' Compute confidence intervals for the data based on some changes. The current
#'  version is tuned to mean changes.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param changes Numeric vector for detected change points.
#' @param K Function for the Kernel. Default is bartlett_kernel.
#' @param h Numeric for bandwidth in computation of long run variance. The default
#'  is \eqn{2N^{1/5}}.
#' @param weighting Weighting for the interval computation, value in \[0,1\].
#'  Default is 0.5.
#' @param M Numeric for the number of Brownian motion simulations in computation
#'  of the confidence interval. Default is 1000.
#' @param alpha Numeric for the significance level, in \[0,1\]. Default is 0.1.
#' @param method String to indicate the method for computing the percentiles used
#'  in the confidence intervals. The options are 'distribution' and 'simulation'.
#'  Default is 'distribution'.
#'
#' @references Horvath, L., & Rice, G. (2024). Change Point Analysis for Time
#'  Series (First edition.). Springer.
#'
#' @references Aue, A., Rice, G., & Sonmez, O. (2018). Detecting and dating structural
#'  breaks in functional data without dimension reduction. Journal of the Royal
#'  Statistical Society. Series B, Statistical Methodology, 80(3), 509-529.
#'
#' @return Data.frame with the first column for the changes, second for the lower
#'  bounds of confidence intervals, and the third for the upper bounds of
#'  confidence intervals.
#' @export
#'
#' @examples
#' X <- cbind(generate_brownian_motion(100,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+1000)
#' confidence_interval(X,changes = 100)
#' confidence_interval(X,changes=100,method = 'simulation')
#'
#' X <- cbind(generate_brownian_motion(100,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.5)
#' confidence_interval(X,100,alpha = 0.1)
#' confidence_interval(X,changes=100,alpha = 0.1,method = 'simulation')
#'
#' X <- generate_brownian_motion(200,v=seq(0,1,0.05))
#' confidence_interval(X,100)
#' confidence_interval(X,100,method = 'simulation')
#'
#' X <- cbind(generate_brownian_motion(200,v=seq(0,1,0.05))$data,
#'            generate_brownian_motion(100,v=seq(0,1,0.05))$data+0.1,
#'            generate_brownian_motion(150,v=seq(0,1,0.05))$data-0.05)
#' confidence_interval(X,c(200,300))
#'
#' confidence_interval(X = electricity, changes = c(64, 120),alpha = 0.1)
confidence_interval <- function(X, changes, K=bartlett_kernel,
                                h=2*ncol(X)^(1/5), weighting=0.5, M=5000,
                                alpha=0.1, method='distribution'){
  ## Verify Inputs
  method <- .verify_input(method, c('distribution','simulation'))
  if(weighting>1 || weighting<0){
    stop('The parameter weighting must be in [0,1].', call. = FALSE)
  }
  if(alpha>1 || alpha<0){
    stop('The parameter alpha must be in [0,1].', call. = FALSE)
  }
  changes <- changes[changes >0]
  changes <- changes[changes <ncol(X)]
  if(length(changes)==0){
    stop('Must specify at least one change in 1, 2, ..., N in parameter changes.', call. = FALSE)
  }
  if(h<0) stop('The parameter h cannot be negative.',call. = FALSE)
  changes <- changes[order(changes)]
  X <- dfts(X)

  changes_ext <- c(0,changes,ncol(X$data))
  r <- nrow(X$data)
  results <- data.frame('change' = rep(NA,length(changes)),
                        'lower' =  NA, 'upper' = NA)

  for(i in 1:length(changes)){
    CP <- changes[i]
    CP_simple <- CP-changes_ext[i]
    idx1 <- (changes_ext[i]+1):CP
    idx2 <- (CP+1):changes_ext[i+2]
    N <- length( c(idx1,idx2) )

    X_k1 <- rowMeans(X$data[,idx1])
    X_k2 <- rowMeans(X$data[,idx2])

    e <- X_k2 - X_k1
    e_l2norm <- dot_integrate(e^2, r = X$fparam)
    eps_hat <- matrix(ncol=N, nrow=r)
    eps_hat[,1:CP_simple] <- X$data[,idx1] - X_k1
    eps_hat[,(CP_simple+1):N] <- X$data[,idx2] - X_k2

    g <- rep(NA,N)
    for(j in 1:N){
      g[j] <- dot_integrate(
        eps_hat[,j] * e / sqrt(e_l2norm),
        r = X$fparam)
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
          m[t_idx] <- (1 - weighting)*(1 - theta) + weighting*theta
        }else if(double_t[t_idx] > 0){
          m[t_idx] <- (1 - weighting)*theta + weighting*(1 - theta)
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

      thetas0 <- stats::quantile(xi, probs = c(alpha/2, 1-alpha/2))

    } else if(method=='distribution'){
      # Pg 40 (pdf 52), information
      .h_func <- function(u,x,y){
        2*x * (x+2*y) * (1- stats::pnorm( (x+2*y)*sqrt(u) )) *
          exp(2*y * (x+y) * u) - 2*x^2 * (1-stats::pnorm( x*sqrt(u) ))
      }
      .f_distribution <- function(t, weighting, theta){
        ifelse(t<=0,
          suppressWarnings(.h_func(-t, (1-weighting)*(1-theta)+weighting*theta, (1-weighting)*theta + weighting*(1-theta))),
          suppressWarnings(.h_func(t, (1-weighting)*theta+weighting*(1-theta), (1-weighting)*(1-theta) + weighting*theta))
        )
      }

      thetas0 <- rep(NA,2)
      f_CDF <- function(up, weighting, theta, alpha){
        abs(stats::integrate(.f_distribution,-20,up, weighting=weighting, theta=theta)$value - alpha)
      }
      # thetas0[1] <- -10.63
      #   # dot_integrate(.f_distribution(seq(-25,-10.63,0.01),weighting,theta),
      #   #                      seq(-25,-10.63,0.01))
      # thetas0[2] <- 11.4
      #   # dot_integrate(.f_distribution(seq(-25,11.4,0.01),weighting,theta),
      #   #                      seq(-25,11.4,0.01))

      tmp <- stats::optimize(f_CDF, interval = c(-20,20),
                      weighting=weighting, theta=theta, alpha=alpha/2)$minimum
      thetas0[1] <- tmp#.f_distribution(tmp,weighting,theta)
      tmp <- stats::optimize(f_CDF, interval = c(-20,20),
                      weighting=weighting, theta=theta, alpha=1-alpha/2)$minimum
      thetas0[2] <- tmp#.f_distribution(tmp,weighting,theta)

    }else{
      stop('Method incorrectly specified',call. = FALSE)
    }

    results[i,] <- c(CP,
                     as.numeric( CP + tau2*thetas0[1] / e_l2norm ),
                     as.numeric( CP + tau2*thetas0[2] / e_l2norm ) )

  }

  results
}

