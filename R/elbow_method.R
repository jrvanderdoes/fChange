#' Elbow Method
#'
#' Method to determine the number of change points using the elbow method. Note,
#'     cascading change points are not considered to allow for every possible
#'     number of change points to be selectable.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param method Method of detecting changes. See [fchange()].
#' @param W Basis to measure the space with when using a characteristic-based
#'  approach.
#' @param trim_function Function with data input. Used to trim data.
#' @param max_changes Numeric for the maximum number of change points to
#'  investigate. Default is the minimum between the data length and 20 changes
#' @param errors String of 'L2' or 'Trace' indicating the error function to use
#' @param K Kernel using for long run covariance
#' @param d Integer for which eigenvalue to compare
#' @param h Lag for long run covariance
#' @param weighting Weights for weighted statistics
#' @param recommendation_change_points Integer for the number of nodes to look ahead at and
#'  confirm a flattening. Default is 2.
#' @param recommendation_improvement Numeric to quantify the gain required of each additional
#'  change. Default is 0.1.
#' @param TVE Total variance explained for projection methods
#'
#' @return List with three elements.
#'  \itemize{
#'    \item **CPInfo**: Data.frame of candidate changes, ordered by impact. The
#'      columns are: CP, Var, and Percent. CP indicates change location, Var is
#'      total variance, Percent is the percentage of variance removed compared
#'      to the previous step. The first row has CP=NA, indicating the no change
#'      scenario.
#'    \item **VarPlot**: A ggplot object showing the variance for increasing
#'      number of changes.
#'    \item **PerPlot**: A ggplot object showing the percent of variance removed
#'      compared to the previous step.
#'  }
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item results <- .elbow_method( X = electricity$data[, 1:50],
#'                                    trim_function = function(...) { 10 },
#'                                    max_changes = 2)
#'  }
.elbow_method <- function(X, method, W=NULL,
                         trim_function = function(X) { 0 },
                         max_changes = min(ncol(X),20),
                         errors = "L2",
                         K=bartlett_kernel,
                         d=NULL, h=0, weighting=0,
                         recommendation_change_points = 2,
                         recommendation_improvement = 0.15,
                         TVE=0.95) {

  ## Setup
  X <- dfts(X)
  max_changes <- min(max_changes, ncol(X))

  n <- ncol(X$data)
  return_data <- data.frame(
    "CP" = NA,
    "Var" = NA
  )

  ## Trim & stopping criteria
  trim_amt <- trim_function(X$data)#, ...)
  if (1 + trim_amt >= n - trim_amt) return()

  ## Run First CP Detection
  if( method=='characteristic' && is.null(W)){
    stop('W must not be null when using this method',call. = FALSE)
  }

  stats_tmp <- .elbow_statistics(X$data[, (1 + trim_amt):(n - trim_amt), drop=FALSE],
                                 method=method,
                                 W=W, d=d, h=h, K=K, weighting=weighting, TVE=TVE)
  test_stats <- c(rep(NA, trim_amt), stats_tmp, rep(NA, trim_amt))
  location <- which.max(stats_tmp) + trim_amt

  ## Total Var Remaining
  return_data[1, ] <- c(location,
                        .compute_total_var(X = X, changes = location,
                                           errors = errors, K=K,
                                           W=W))

  ## Iteratively Search
  while (nrow(return_data) < max_changes) {
    ## Setup
    changes <- c(0, return_data$CP, n)
    changes <- changes[order(changes)]
    prev_CP <- return_data[nrow(return_data), "CP"]
    prev_CP_loc <- which(changes == prev_CP)

    ## We only need to recompute for the interval changed by last CP!
    # Before
    beforePrevCP <- (changes[prev_CP_loc - 1] + 1):(changes[prev_CP_loc])
    trim_amt <- trim_function(X$data[, beforePrevCP])#, ...)
    test_stats[beforePrevCP] <- 0

    if (1 + trim_amt < length(beforePrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(beforePrevCP) - trim_amt)
      stats_tmp <- .elbow_statistics(X$data[, beforePrevCP[fill_idx], drop=FALSE],
                                     method=method,
                                     W=W, d=d, h=h, K=K, weighting=weighting,
                                     TVE=TVE)

      # Set CP to 0 if include (i.e. not trimmed)
      if(max(beforePrevCP[fill_idx])==prev_CP) stats_tmp[length(stats_tmp)] <- 0

      test_stats[beforePrevCP[fill_idx]] <- stats_tmp
    }

    # After
    afterPrevCP <- (changes[prev_CP_loc] + 1):changes[prev_CP_loc + 1]
    trim_amt <- trim_function(X$data[, afterPrevCP])#, ...)
    test_stats[afterPrevCP] <- 0

    if (1 + trim_amt < length(afterPrevCP) - trim_amt) {
      fill_idx <- (1 + trim_amt):(length(afterPrevCP) - trim_amt)
      stats_tmp <- .elbow_statistics(X$data[, afterPrevCP[fill_idx], drop=FALSE],
                                 method=method,
                                 W=W, d=d, h=h, K=K, weighting=weighting, TVE=TVE)

      # Set next CP to 0 if include (i.e. not trimmed)
      if(max(afterPrevCP[fill_idx])==changes[prev_CP_loc + 1]) stats_tmp[length(stats_tmp)] <- 0

      test_stats[afterPrevCP[fill_idx]] <- stats_tmp
    }

    #### Split based on 0s (changes and trimmings)
    is.sep <- test_stats==0
    data_segments <- split(test_stats[!is.sep], cumsum(is.sep)[!is.sep])


    ## Get total variance for each potential CP
    return_data_tmp <- data.frame("CP" = rep(NA,length(data_segments)), "Var" = NA)
    for (k in 1:length(data_segments)) {
      ## Find max in section
      value_max <- max(data_segments[[k]])
      section_max <- min(which(test_stats == value_max))
      changes_proposed <- c(section_max, return_data$CP)

      return_data_tmp[k, ] <- c(
        section_max,
        .compute_total_var(X = X, changes = changes_proposed, errors = errors, W=W, K=K)
      )
    }

    ## With max test-statistic on each section, take one leading to min variance
    return_data[nrow(return_data)+1, ] <- return_data_tmp[which.min(return_data_tmp$Var), ]

    if(sum(test_stats,na.rm = TRUE)==0) break
  }


  ## Add No Change option and compute percent change
  return_data <- rbind(
    c(NA, .compute_total_var(X=X, changes=c(), errors=errors, W=W, K=K)),
    return_data
  )
  return_data$Percent <- 1 - return_data$Var / max(return_data$Var)
  return_data$Improvement <- c(NA,1-return_data$Var[-1] /return_data$Var[-nrow(return_data)])

  ## Recommended Number of changes
  #   Look ahead at the next recommendation_change_points changes
  #   Determine if the reduction in variance is less than requested
  #   If so, select. Otherwise look at next point
  data_change<- data.frame(matrix(NA,nrow=nrow(return_data),
                                  ncol = recommendation_change_points))
  for(i in 1:recommendation_change_points){
    data_change[,i] <-
      c(return_data$Improvement[-c(1:i)]<recommendation_improvement, rep(TRUE,i))
  }
  recommended_cp <- which.max(apply(data_change, MARGIN = 1, prod)) - 1
  if(recommended_cp==0){
    recommended_changes <- NA
  } else{
    recommended_changes <- return_data$CP[1:(recommended_cp+1)]
    recommended_changes <- recommended_changes[!is.na(recommended_changes)]
  }


  ## Define vars to remove notes
  CP <- Var <- Percent <- Improvement <- NULL

  var_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Variability") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  per_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent),size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Percent),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Percent of Total Variability Explained") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  ## TODO:: CHECK THIS OUT
  gain_plot <-
    ggplot2::ggplot(return_data[-1,]) +
    ggplot2::geom_point(ggplot2::aes(x = 1:length(CP), y = Improvement), size=4) +
    ggplot2::geom_line(ggplot2::aes(x = 1:length(CP), y = Improvement),linewidth=2) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Variability Improvement from Previous") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))


  recommend_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_line(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),
                       linewidth=2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=recommended_cp),
                        linetype='dotted', col='red',linewidth=2) +
    ggplot2::geom_point(ggplot2::aes(x = 0:(length(CP) - 1), y = Var),
                        size=4) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Total Variance") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                   axis.title = ggplot2::element_text(size=22))

  list("information" = return_data,
       "plots"=list(
         "variability" = var_plot,
         "explained" = per_plot,
         "improvement" = gain_plot
         ),
       "suggestion" = list(
         'plot'=recommend_plot,
         'changes'=recommended_changes)
       )
}



#' Compute Test Statistics for Elbow Plot
#'
#' @inheritParams .elbow_method
#'
#' @returns Statistic based on the method given
#'
#' @noRd
#' @keywords internal
.elbow_statistics <- function(X, method, W=NULL, d=NULL, h=NULL, K=NULL,
                              weighting=0, TVE=0.95){

  n <- ncol(X)
  r <- nrow(X)

  if(method=='characteristic'){

    # Compute Zn
    fhat_vals <- as.matrix(.fhat_all(X, W))
    Zn <- sqrt(n) * (fhat_vals - (seq_len(n)/n) %o% fhat_vals[nrow(fhat_vals),])

    # ns <- .select_n(1:n, J)
    # Zn <- Zn[ns,,drop=FALSE]

    # Integrate out W
    stats <- dot_integrate_col(t(abs(Zn)^2))
  } else if(method=='mean'){
    stats <- rep(0, n)
    # CUSUM
    for (k in 1:n) {
      stats[k] <- sum((rowSums(X[, 1:k,drop=FALSE]) -
                       (k / n) * rowSums(X))^2)
    }
    stats <- stats / n
  }  else if(method=='robustmean'){
    # U_N(x) stacked
    hC_Obs <- make_hC_Obs(X)

    # find max_k U_{n,k} (missing constants)
    stats <- c(fill_T(hC_Obs, n, r)[,1],0)

  } else if(method=='eigenjoint'){

    Cov_op <- .partial_cov(X, 1)

    eig2 <- array(dim = c(r,r,r))
    for(j in 1:r){
      eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
    }
    thetas <- matrix(nrow=n,ncol=r)
    X2 <- array(dim = c(r,r,n))
    for(i in 1:n){
      X2[,,i] <- X[,i] %*% t(X[,i]) - Cov_op$coef_matrix

      for(j in 1:r){
        thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
      }
    }

    Sigma_d <- long_run_covariance(X = t(thetas[,1:d]), h = h, K = K)
    # Moore Penrose solve if non-invertable
    Sigma_d_inv <- tryCatch({
      solve(Sigma_d)
    }, error = function(e){
      eigs <- eigen( Sigma_d )
      eigs$vectors <- eigs$vectors * sqrt(r)
      eigs$values <- eigs$values / r

      K_eigs <- which.min(cumsum(eigs$values)/sum(eigs$values) < 0.95)

      tmp_inv <- matrix(0, nrow=nrow(Sigma_d), ncol=ncol(Sigma_d))
      for(i in 1:K_eigs){
        tmp_inv <- tmp_inv +
          ( 1 / eigs$values[i] ) * (eigs$vectors[,i] %o% eigs$vectors[,i] )
      }

      tmp_inv
    })


    stats <-  rep(0,n)
    for(k in 1:n){
      lam <- .partial_cov(X, k/n)$eigen_val[1:d]
      kapa <- lam - (k/n) * Cov_op$eigen_val[1:d]
      stats[k] <- n * t(kapa) %*% Sigma_d_inv %*% kapa
    }
  } else if(method=='eigensingle'){

    Cov_op <- .partial_cov(X, 1)

    eig2 <- array(dim = c(r,r,r))
    for(j in 1:r){
      eig2[,,j] <- Cov_op$eigen_fun[,j] %*% t(Cov_op$eigen_fun[,j])
    }
    thetas <- matrix(nrow=n,ncol=r)
    X2 <- array(dim = c(r,r,n))
    for(i in 1:n){
      X2[,,i] <- X[,i] %*% t(X[,i]) - Cov_op$coef_matrix

      for(j in 1:r){
        thetas[i,j] <- sum(diag( t(eig2[,,j]) %*% X2[,,i] ))
      }
    }

    Sigma_d <- long_run_covariance(X = t(thetas[,d]), h = h, K = K)

    stats <- rep(0,n)
    for (k in 1:n){
      lam <- .partial_cov(X, k/n)$eigen_val[d]
      stats[k] <- ( n*( lam - (k/n)*Cov_op$eigen_val[d] )^2 ) / Sigma_d
    }

  }else if(method=='trace'){

    Cov_op <- .partial_cov(X, 1)
    lambda <- Cov_op$eigen_val
    T_1 <- sum(lambda)
    Xi <- sapply(1:n, function(i,X1){ X[,i] %*% X[,i] }, X1=X)
    sigma_sq <- tryCatch({
      suppressWarnings(sandwich::lrvar(Xi, prewhite = FALSE))
    },error=function(e){
      NA
    })
    if(is.na(sigma_sq)) return(rep(0, n))

    sigma <- sqrt(sigma_sq)

    stats <- rep(0,n)
    for (k in 1:n) {
      T_x <- sum(Xi[1:k])/n
      stats[k] <- (1/sigma) * abs(T_x - k/n * T_1)
    }
  }else if(method=='covariance'){

    xdm <- center(dfts(X))$data

    uind <- seq(0, 1, length = n + 1)[2:(n + 1)]
    zn2 <- list()
    stats <- c(rep(0, n))
    for (i in 1:(n - 1)) {
      Zn_stat <- .covariance_statistic_cp(xdm, uind[i])

      stats[i] <- (n / (i * (n - i)))^weighting *
        dot_integrate(dot_integrate_col( Zn_stat^2))
    }

  }else if(method=='projmean'){

    pca_X <- pca(dfts(X), TVE=TVE)
    d <- length(pca_X$sdev)

    eta.hat <- as.matrix(pca_X$x)

    ## Test Statistic
    kappa <-
      sapply(1:n,function(k, eta.hat, n){
        colSums(eta.hat[1:k, , drop=FALSE]) - k/n * colSums(eta.hat)
      },eta.hat=eta.hat, n=n)
    # If same names, only 1 pc and sapply flips it incorrectly
    if(length(unique(names(kappa)))==1) kappa <- t(kappa)

    if(length(pca_X$sdev)>1){
      stats <-  diag( 1/n * ( t(kappa) %*% diag(1/pca_X$sdev^2) %*% kappa ) )
    }else{
      stats <-  diag( 1/n * ( t(kappa) %*% (1/pca_X$sdev^2) %*% kappa ) )
    }
  }else if(method=='projdistribution'){

    X_pca <- pca(dfts(X), TVE=TVE)
    n <- ncol(X)
    d <- length(X_pca$sdev)

    kappa <- matrix(nrow=d, ncol=n-1)

    # Params for loop
    t_vals <- seq(0, 1, length.out=n )
    ws <- .w(t_vals)
    ks <- 1:(n-1)
    weights <- ((ks * (n - ks)) / n^2)^weighting * ((ks * (n - ks)) / n)

    for(i in 1:d){
      Y <- X_pca$x[,i]

      dat <- exp(complex(real = 0, imaginary = 1) * t_vals %*% t(Y) )

      ## This commented code is a slow version of below
      # internals <- sapply(ks, function(k, dat, ws){
      #   abs(rowMeans(dat[,1:k,drop=FALSE]) - rowMeans(dat[,(k+1):n,drop=FALSE]) )^2 * ws
      # },dat=dat,ws=ws)

      cmean <- t(apply(dat,MARGIN = 1,cumsum)) /
        matrix(1:ncol(dat),nrow=length(t_vals),ncol=n,byrow = TRUE)
      cmean1 <- t(apply(dat[,n:2,drop=FALSE],MARGIN = 1,cumsum))[,(n-1):1] /
        matrix((n-1):1,nrow=length(t_vals),ncol=n-1,byrow = TRUE)
      internals <- abs(cmean[,-n] - cmean1)^2 * ws

      kappa[i,] <- weights * dot_integrate_col(internals)
    }
    if(d>1){
      stats <-  diag( 1/n * ( t(kappa) %*% diag(1/X_pca$sdev^2) %*% kappa ) )
    }else{
      stats <-  diag( 1/n * ( t(kappa) %*% (1/X_pca$sdev^2) %*% kappa ) )
    }

    stats <- c(stats,0)
  }

  # Make last value 0
  stats[length(stats)] <- 0

  stats
}


#' Compute Total Variance
#'
#' This (internal) function computes the total variance in the data with given
#'  changes.
#'
#' @inheritParams .elbow_method
#'
#' @return Numeric indicating the variance between all subsegments
#'
#' @noRd
#' @keywords internal
.compute_total_var <- function(X, changes, errors = "L2", K=bartlett_kernel,
                               W=generate_brownian_motion(100, v=X$fparam )$data) {
  # Get Information
  if (tolower(errors) == "l2" || tolower(errors) == "trace") {

    X_std <- center(X, changes=changes)$data
    covMatrix <- long_run_covariance(X_std, h=0, K=K)

  }
  # else if (tolower(errors) == "ce") {
  #
  #   CE <- data.frame(matrix(ncol = ncol(X), nrow = ncol(W)))
  #
  #   ## Compute CEs
  #   for (i in 1:ncol(X)) {
  #     CE[, i] <- apply(W, 2,
  #                      FUN = function(v, dat) {
  #                        exp(complex(real = 0, imaginary = 1) * (t(dat) %*% v))
  #                      }, dat = X$data[, i]
  #     )
  #   }
  #
  #   CE_std <- center(dfts(CE),changes=changes)$data
  #   covMatrix <- long_run_covariance(t(CE_std) %*% CE_std, h=0, K=K)
  #
  #
  # }

  # Apply error metric
  if (tolower(errors) == "l2") {
    returnValue <- sqrt(sum(covMatrix^2))
  } else if (tolower(errors) == "trace") {
    returnValue <- sum(diag(covMatrix))
  } else {
    stop("Only L2 and Trace error functions implemented")
  }
  # else if (tolower(errors) == "ce") {
  #   # returnValue <- sum(colSums(abs(CE_std)^2) / ncol(W))
  #   returnValue <- sqrt(sum(abs(covMatrix)^2))
  # }

  returnValue
}



