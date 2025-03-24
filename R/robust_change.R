#' Robust Change Point Detection
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param m Number of interations for resampling
#' @param statistic Test statistic of interest. Either Integrated (\code{Tn}) or
#'  maximized (\code{Mn}). Default is \code{Tn}
#' @param statistic Threshold method of interest. Either the simulation distribution
#'  (\code{simulation}) or resampling (\code{resample}). Default is simulation.
#'
#' @return A list with the following elements
#' \itemize{
#'  \item pvalue: P-value for the test statistic
#'  \item statistic: Test statistic value
#'  \item simulations: Values from the distribution simulations or permuated
#'    simulations
#'  \item bandwidth: Bandwidth for test
#' }
#'
#' @noRd
#' @keywords internal
#'
#' @references Wegner, L., Wendler, M. Robust change-point detection for
#'  functional time series based on U-statistics and dependent wild bootstrap.
#'  Stat Papers (2024).
#'
#' @examples
#' #result <- .change_robust(dfts(electricity$data[,1:100]),m=10)
.change_robust <- function(X, m, statistic=c('Tn','Mn'),
                          threshold=c('simulation','resample'),
                          bandwidth = NA){
  # Setup variables
  statistic <- .verify_input(statistic, c('Tn','Mn'))
  threshold <- .verify_input(threshold, c('simulation','resample'))
  X <- dfts(X)

  # Create Proper Data and Bandwidths
  Obs <- X$data

  ## Get Information - Cov Matrix
  Obs_tilde_h <- make_Obs_tilde_h(Obs)
  if(is.na(bandwidth)){
    bw_h <- adaptive_bw(Obs_tilde_h)
  }else{
    bw_h <- bandwidth
  }

  # Calculate quadratic spectral cov. matrices
  ker_h <- kernel(bw_h, ncol(Obs))
  A_h <- toeplitz(ker_h)
  Re_A_h <- getRealSQM(A_h)

  ##########
  ## Compute Test Statistic And Threshold
  ##########
  # U_N(x) stacked
  hC_Obs <- make_hC_Obs(Obs)

  # find max_k U_{n,k} (missing constants)
  Te <- fill_T(hC_Obs, ncol(Obs), nrow(Obs))

  if(statistic=='Tn'){
    stat <- dot_integrate(Te^2) / ncol(Obs)^(3)
  }else if(statistic=='Mn'){
    # T_max <- apply(Te,MARGIN = 2, max)
    stat <- max(Te) / ncol(Obs)^(3/2)
  }
  location <- which.max(Te)

  ##########################
  # # MUCH SLOWER, JUST FOR UNDERSTANDING
  # h_func <- function(x,y){ (x-y) / sqrt( sum((x-y)^2) ) }

  # ## Examine hC_Obs, T_Fill, and This
  # N <- ncol(Obs)
  # UN <- matrix(NA, nrow = nrow(Obs), ncol = N-1)
  # for(k in 1:(N-1)){
  #   UN_tmp <- matrix(0, nrow=3, ncol=N-k )
  #   for(i in 1:k){
  #     for(j in (k+1):N){
  #       UN_tmp[,j-k] <- UN_tmp[,j-k] + h_func(Obs[,i], Obs[,j])
  #     }
  #   }
  #
  #   UN[,k] <- rowSums(UN_tmp)
  # }
  # stat = max( sqrt(colSums(UN^2)) )

  # # Obs_tilde_h
  # h1_X <- matrix(0,nrow = nrow(Obs), ncol = ncol(Obs))
  # for(i in 1:ncol(Obs)){
  #   for(j in 1:ncol(Obs)){
  #     if(j==i)
  #       next
  #     h1_X[,i] <- h1_X[,i] + h_func(Obs[,i], Obs[,j])
  #   }
  # }
  # tmp <- 1/(ncol(Obs)-1)* h1_X


  # ######################
  # ### simulation (TN)
  # dist_values <- dist_values1 <- dist_values2 <- dist_values3 <-
  #   dist_values4 <- dist_values5 <-
  #   dist_values6 <- dist_values7 <-
  #   dist_values8 <- dist_values9 <-
  #   dist_values10 <- dist_values11 <- rep(NA, m)
  # b_res <- 250
  # for(i in 1:m){
  #   BridgeLam <- BridgeLam1 <- BridgeLam2 <- BridgeLam3 <-
  #     BridgeLam4 <- BridgeLam5 <-
  #     rep(NA,length.out=length(lambda))
  #   for(j in 1:length(lambda)){
  #     BridgeLam[j] <- lambda[j] *
  #       dot_integrate(sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #     BridgeLam1[j] <- lambda[j] *
  #       dot_integrate(sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #     BridgeLam2[j] <- lambda1[j] *
  #       dot_integrate(sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #     BridgeLam3[j] <- lambda1[j] *
  #       dot_integrate(sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #     BridgeLam4[j] <- lambda2[j] *
  #       dot_integrate(sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #     BridgeLam5[j] <- lambda2[j] *
  #       dot_integrate(sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2)
  #   }
  #   dist_values[i] <- sum(BridgeLam)
  #   dist_values1[i] <- sum(BridgeLam1)
  #   dist_values2[i] <- sum(BridgeLam2)
  #   dist_values3[i] <- sum(BridgeLam3)
  #   dist_values4[i] <- sum(BridgeLam4)
  #   dist_values5[i] <- sum(BridgeLam5)
  #
  #   dist_values6[i] <- sqrt(sum(BridgeLam^2))
  #   dist_values7[i] <- sqrt(sum(BridgeLam1^2))
  #   dist_values8[i] <- sqrt(sum(BridgeLam2^2))
  #   dist_values9[i] <- sqrt(sum(BridgeLam3^2))
  #   dist_values10[i] <- sqrt(sum(BridgeLam4^2))
  #   dist_values11[i] <- sqrt(sum(BridgeLam5^2))
  # }
  #
  # mean(stat <= dist_values)
  # mean(stat <= dist_values1)
  # mean(stat <= dist_values2)
  # mean(stat <= dist_values3) #
  # mean(stat <= dist_values4) #
  # mean(stat <= dist_values5) #
  #
  # mean(stat <= dist_values6)
  # mean(stat <= dist_values7)
  # mean(stat <= dist_values8)
  # mean(stat <= dist_values9)
  # mean(stat <= dist_values10) #
  # mean(stat <= dist_values11) #

  #######################
  ### simulation (MN)
  # dist_values <- dist_values1 <- dist_values2 <- dist_values3 <-
  #   dist_values4 <- dist_values5 <-
  #   dist_values6 <- dist_values7 <-
  #   dist_values8 <- dist_values9 <-
  #   dist_values10 <- dist_values11 <- rep(NA, m)
  # b_res <- 250
  # for(i in 1:m){
  #   BridgeLam <- BridgeLam1 <- BridgeLam2 <- BridgeLam3 <-
  #     BridgeLam4 <- BridgeLam5 <-
  #     matrix(NA,nrow=length(lambda1), ncol=b_res)
  #   for(j in 1:length(lambda1)){
  #     # BridgeLam[j,] <- lambda[j] *
  #     #   sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2
  #     # BridgeLam1[j,] <- lambda[j] *
  #     #   sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2
  #     BridgeLam2[j,] <- lambda1[j] *
  #       sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2
  #     # BridgeLam3[j,] <- lambda1[j] *
  #     #   sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2
  #     # BridgeLam4[j,] <- lambda2[j] *
  #     #   sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2
  #     # BridgeLam5[j,] <- lambda2[j] *
  #     #   sde::BM(x = 0, t0 = 0, T = 1, N = b_res-1)^2
  #   }
  #   # dist_values[i] <- max(sqrt(colSums(BridgeLam)))
  #   # dist_values1[i] <- max(sqrt(colSums(BridgeLam1)))
  #   dist_values2[i] <- max(sqrt(colSums(BridgeLam2)))
  #   # dist_values3[i] <- max(sqrt(colSums(BridgeLam3)))
  #   # dist_values4[i] <- max(sqrt(colSums(BridgeLam4)))
  #   # dist_values5[i] <- max(sqrt(colSums(BridgeLam5)))
  #   #
  #   # dist_values6[i] <- max(sqrt(sqrt(colSums(BridgeLam^2))))
  #   # dist_values7[i] <- max(sqrt(sqrt(colSums(BridgeLam1^2))))
  #   # dist_values8[i] <- max(sqrt(sqrt(colSums(BridgeLam2^2))))
  #   # dist_values9[i] <- max(sqrt(sqrt(colSums(BridgeLam3^2))))
  #   # dist_values10[i] <- max(sqrt(sqrt(colSums(BridgeLam4^2))))
  #   # dist_values11[i] <- max(sqrt(sqrt(colSums(BridgeLam5^2))))
  #
  #   # dist_values6[i] <- max(sqrt((colSums(BridgeLam^2))))
  #   # dist_values7[i] <- max(sqrt((colSums(BridgeLam1^2))))
  #   # dist_values8[i] <- max(sqrt((colSums(BridgeLam2^2))))
  #   # dist_values9[i] <- max(sqrt((colSums(BridgeLam3^2))))
  #   # dist_values10[i] <- max(sqrt((colSums(BridgeLam4^2))))
  #   # dist_values11[i] <- max(sqrt((colSums(BridgeLam5^2))))
  # }
  #
  # return(
  #   c(
  #     # mean(stat <= dist_values),
  #     # mean(stat <= dist_values1),
  #     mean(stat <= dist_values2)#,
  #     # mean(stat <= dist_values3),
  #     # mean(stat <= dist_values4),
  #     # mean(stat <= dist_values5),
  #     #
  #     # mean(stat <= dist_values6),
  #     # mean(stat <= dist_values7),
  #     # mean(stat <= dist_values8),
  #     # mean(stat <= dist_values9),
  #     # mean(stat <= dist_values10),
  #     # mean(stat <= dist_values11)
  #   )
  # )
  #
  # mean(stat <= dist_values)
  # mean(stat <= dist_values1)
  # mean(stat <= dist_values2)
  # mean(stat <= dist_values3)
  # mean(stat <= dist_values4)
  # mean(stat <= dist_values5)
  #
  # mean(stat <= dist_values6)
  # mean(stat <= dist_values7)
  # mean(stat <= dist_values8)
  # mean(stat <= dist_values9)
  # mean(stat <= dist_values10)
  # mean(stat <= dist_values11)

  #########################

  if(threshold=='simulation'){
    ## TODO: Set this to use same kernel
    cov_mat <- long_run_covariance(X = Obs_tilde_h, h = bw_h, K = bartlett_kernel)
    lambda <- eigen(cov_mat)$values

    sims <- rep(NA, m)
    b_res <- 250

    if(statistic=='Tn'){
      for(i in 1:m){
        BridgeLam <- lambda *
          dot_integrate_col(
            generate_brownian_bridge(length(lambda),v=seq(0,1,length.out=b_res))$data^2
          )
        # BridgeLam <- rep(NA,length(lambda))
        # for(j in 1:length(lambda)){
        #   BridgeLam[j] <- lambda[j] *
        #     dot_integrate( sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2 )
        # }

        sims[i] <- sum(BridgeLam)
      }

    }else if(statistic=='Mn'){
      for(i in 1:m){
        # BridgeLam <- matrix(NA,nrow=length(lambda), ncol=b_res)
        # for(j in 1:length(lambda)){
        #   BridgeLam[j,] <- lambda[j] *
        #     sde::BBridge(x = 0, y = 0, t0 = 0, T = 1, N = b_res-1)^2
        # }
        BridgeLam <-
          generate_brownian_bridge(length(lambda),v=seq(0,1,length.out=b_res))$data^2 %*%
          diag(lambda)
        sims[i] <- max(sqrt(rowSums(BridgeLam)))
      }
    }

  }else if(threshold=='resample'){
    # MULTIPLIER BOOTSTRAP
    #   m x n
    U_hC <- fill_U(A_h=Re_A_h,
                   # A_C=Re_A_C,
                   norm=stats::rnorm(m*ncol(Obs)),
                   hC_Obs=hC_Obs,
                   m=m, n=ncol(Obs), d=nrow(Obs) )

    # Get Sup at each Sim
    if(statistic=='Tn'){
      sims <- dot_integrate_col(t(U_hC)^2) / ncol(Obs)^(3)

    }else if(statistic=='Mn'){
      # find max_k U_{n,k}^(t) for t=1,...,m
      sims <- 1/ncol(Obs)^(3/2) * find_rowmax(U_hC)
    }

  }


  list("pvalue" = sum(stat <= sims)/m,
       "location" = location)#,
       # 'statistic' = stat,
       # 'simulations' = sims,
       # 'extra' = list("bandwidth"=bw_h)
       #)
}
