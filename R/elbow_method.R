
#' Elbow Method
#'
#' Method to determine the number of change points using the elbow method. Note,
#'     cascading change points are not considered to allow for every possible
#'     number of change points to be selectable.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param test_statistic_function Function with the first argument being data
#'     and the second argument argument for candidate change points.
#'     Additional arguments passed in via ... . Return a single numeric value.
#' @param cutoff_function Function with first argument being data and the second
#'     argument being alpha. No other arguments given. Return single numeric
#'     value.
#' @param trim_function Function taking data as an argument and returning a
#'     numeric value indicating how much should be trimmed on each end
#' @param alpha Numeric value in [0,1] indicating the significance for
#'     cutoff_function.
#' @param errorType String of 'L2' or 'Tr' indicating the error function to use
#' @param ... Additional parameters to pass into the respective functions
#'
#' @return list with element 1 the data frame with the change point location and
#'     element a ggplot of variance as a function of CPs
#' @export
#'
#' @examples
#' data_KL <- generate_data_fd(ns = c(100,50,100),
#'                   eigsList = list(c(3,2,1,0.5),
#'                                   c(3,2,1,0.5),
#'                                   c(3,2,1,0.5)),
#'                   basesList = list(
#'                      fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'                   meansList = c(-1,1,-1),
#'                   distsArray = c('Normal','Normal','Normal'),
#'                   evals = seq(0,1,0.05),
#'                   kappasArray = c(0, 0, 0))
#'
#' results <- elbow_method(data=data_KL,
#'                   test_statistic_function=compute_Mn, fn=compute_Mn,
#'                   cutoff_function = generalized_resampling,
#'                   trim_function = trim_function,
#'                   alpha = 0.05, errorType = 'L2',
#'                   M=1000)
#' print(results[[2]])
#' plot_fd(data_KL, CPs=results[[1]]$CP[2:3])
elbow_method <- function(data,
                         test_statistic_function =NULL,
                         changepoint_function =NULL,
                         errorType = 'L2', ... ){

  if(!is.null(test_statistic_function)){
    elbow_result <- .elbow_method_stat(data,
                                       test_statistic_function=test_statistic_function,
                                       errorType=errorType, ... )
  }else if(!is.null(changepoint_function)){
    elbow_result <- .elbow_method_change(data,
                                         changepoint_function=changepoint_function,
                                         errorType=errorType, ... )
  }

}

.elbow_method_stat <- function(data, test_statistic_function,
                         cutoff_function, trim_function,
                         alpha=0.05, errorType = 'L2', ... ){
  # Setup
  n <- ncol(data)
  return_data <- data.frame('CP'=NA,
                            'Var'=NA)

  # Trim & stopping criteria
  trim_amt <- trim_function(data)#, ...)
  nStart <- 1+trim_amt
  nEnd <- ncol(as.data.frame(data))-trim_amt
  if(nStart> nEnd) return()

  # Run First
  test_stat <- rep(NA, n)
  for(k in nStart:nEnd){
    test_stat[k] <- test_statistic_function(data, k, ...)
  }
  return_data[1,] <- c(which.max(test_stat),
                       .compute_total_var(data, which.max(test_stat),
                                          errorType))

  # Iteratively Search
  i <- 1
  while(TRUE){
    # Setup
    i <- i + 1
    CPs <- c(0,return_data$CP, n)
    CPs <- CPs[order(CPs)]
    prev_CP <- return_data[nrow(return_data),'CP']
    prev_CP_loc <- which(CPs==prev_CP)

    # Previous
    prev2currCP <- (CPs[prev_CP_loc-1]+1):(CPs[prev_CP_loc])
    for(k in prev2currCP){
      if(k < min(prev2currCP)+trim_amt ||
         k > max(prev2currCP)-trim_amt){
        test_stat[k] <- NA
      }
      else{
        test_stat[k] <- test_statistic_function(
                                as.data.frame(data[, prev2currCP]),
                                k-min(prev2currCP)+1, ...)
      }
    }
    # Next
    curr2nextCP <- (CPs[prev_CP_loc]+1):CPs[prev_CP_loc+1]
    for(k in curr2nextCP){
      if(k < min(curr2nextCP)+trim_amt ||
         k > max(curr2nextCP)-trim_amt){
        test_stat[k] <- NA
      }
      else{
        test_stat[k] <- test_statistic_function(
                                as.data.frame(data[, curr2nextCP]),
                                k-min(curr2nextCP)+1, ...)
      }
    }

    if(sum(is.na(test_stat))==n)
      break

    ## Get total variance for each potential CP
    data_segments <- .split_on_NA(test_stat)

    return_data_tmp <- data.frame('CP'=NA,'Var'=NA)
    for(k in 1:length(data_segments)){
      # Find max test stat on interval
      value_max <- max(data_segments[[k]])

      # Get CP and total variance with full data
      section_max <- which(test_stat==value_max)
      tmp <- c(section_max, return_data$CP)
      tmp <- tmp[order(tmp)]

      return_data_tmp[k,] <- c(section_max,
                               .compute_total_var(data, tmp, errorType))
    }

    ## With max test-statistic on each section, take one leading to min variance
    return_data[i,] <- return_data_tmp[which.min(return_data_tmp$Var),]


    if(nrow(return_data)==(n-1))
      break
  }

  return_data <- rbind(c(NA, .compute_total_var(data, c(), errorType)),
                       return_data)
  return_data$Percent <- 1-return_data$Var/max(return_data$Var)

  var_plot <-
    ggplot(return_data) +
    geom_point(aes(x=0:(length(CP)-1), y=Var)) +
    geom_line(aes(x=0:(length(CP)-1), y=Var)) +
    xlab("Number of Change Points") +
    ylab("Total Variance") +
    theme_bw()

  per_plot <-
    ggplot(return_data) +
    geom_point(aes(x=0:(length(CP)-1), y=Percent)) +
    geom_line(aes(x=0:(length(CP)-1), y=Percent)) +
    xlab("Number of Change Points") +
    ylab("Percent Explained") +
    theme_bw()

  return(list(return_data, var_plot, per_plot))
}


.elbow_method_change <- function(data, changepoint_function, errorType = 'L2', ... ){
  # Setup
  n <- ncol(data)
  return_data <- data.frame('CP'=NA,
                            'Var'=NA)

  # Trim & stopping criteria
  elbow_res <- changepoint_function(data, ...)
  return_data[1,] <- c(elbow_res,
                       .compute_total_var(data, elbow_res, errorType))

  # Iteratively Search
  i <- 1
  while(TRUE){
    # Setup
    i <- i + 1
    CPs <- c(0,return_data$CP, n)
    CPs <- CPs[order(CPs)]

    # Potential CPs
    return_data_tmp <- data.frame('CP'=rep(NA, length(CPs)-1),'Var'=NA)
    for(j in 1:(length(CPs)-1)){
      if(((CPs[j+1]-1)-CPs[j])<=1)
        next
      tmp <- changepoint_function(as.data.frame(data[,(CPs[j]+1):CPs[j+1]]), ...) +
        CPs[j]
      if(is.na(tmp) || (tmp == CPs[j]) || (tmp == CPs[j+1]))
        next
      tmp1 <- c(CPs,tmp)
      tmp1 <- tmp1[order(tmp1)]

      return_data_tmp[j,] <-
              c(tmp, .compute_total_var(data, tmp1[-c(1,length(tmp1))], errorType))
    }

    # Stop if no CP detected
    return_data_tmp <- na.omit(return_data_tmp)
    if(nrow(return_data_tmp)==0)
      break

    ## With max test-statistic on each section, take one leading to min variance
    return_data[i,] <- return_data_tmp[which.min(return_data_tmp$Var),]

    # Stop if all CPs detected
    if(nrow(return_data)==(n-1))
      break
  }

  return_data <- rbind(c(NA, .compute_total_var(data, c(), errorType)),
                       return_data)
  return_data$Percent <- 1-return_data$Var/max(return_data$Var)

  var_plot <-
    ggplot(return_data) +
    geom_point(aes(x=0:(length(CP)-1), y=Var)) +
    geom_line(aes(x=0:(length(CP)-1), y=Var)) +
    xlab("Number of Change Points") +
    ylab("Total Variance") +
    theme_bw()

  per_plot <-
    ggplot(return_data) +
    geom_point(aes(x=0:(length(CP)-1), y=Percent)) +
    geom_line(aes(x=0:(length(CP)-1), y=Percent)) +
    xlab("Number of Change Points") +
    ylab("Percent Explained") +
    theme_bw()

  return(list(return_data, var_plot, per_plot))
}

.compute_total_var <- function(data, CPs, errorType='L2', M=1000){
  # Setup
  data <- as.data.frame(data)
  data_std <- data
  CPs <- c(0, CPs, ncol(data))

  if(errorType=='L2' || errorType=='Tr'){

    # Standardize Data
    for(i in 2:length(CPs)){
      indices <- (CPs[i-1]+1):CPs[i]
      CP_mean <- rowMeans(as.data.frame(data[,indices])) ##TODO:: Verify
      # Standardize Data
      data_std[,indices] <- data_std[,indices]-CP_mean
    }

    ## Cannot use cov because it removes the mean from individual FDs, thus
    #     making each change point near equally effective
    covMatrix <- matrix(NA,ncol=ncol(data), nrow = ncol(data))
    sampleCoef <- 1/(nrow(data)-1)

    for(i in 1:ncol(data)){
      for(j in i:ncol(data)){

        covMatrix[i,j] <- covMatrix[j,i] <-
          sum(sampleCoef*(data_std[,i] * data_std[,j]))
      }
    }
  } else if(errorType=='CE'){

    W <- as.data.frame(sapply(rep(0,M),sde::BM, N=nrow(data)-1))
    CE <- data.frame(matrix(ncol=ncol(data), nrow=M))

    ## Compute CEs
    for(i in 1:ncol(data)){
      CE[,i] <- apply(W, 2,
                      FUN = function(v, dat){
                        exp(complex(real=0,imaginary = 1) * (t(dat) %*% v))
                      }, dat=data[,i])
    }

    CE_std <- CE

    # Standardize Data
    for(i in 2:length(CPs)){
      indices <- (CPs[i-1]+1):CPs[i]
      CP_mean <- rowMeans(as.data.frame(CE[,indices])) ##TODO:: Verify
      # Standardize Data
      CE_std[,indices] <- CE[,indices]-CP_mean
    }

    # covMatrix <- matrix(NA,ncol=ncol(data), nrow = ncol(data))
    # sampleCoef <- 1/(nrow(data)-1)
    #
    # ## Compute covMatrix
    # for(i in 1:ncol(data)){
    #   for(j in (i):ncol(data)){
    #     covMatrix[i,j] <- covMatrix[j,i] <-
    #       sum(sampleCoef*(CE[,i] %*% CE[,j]))
    #   }
    # }
  }

  # Apply error metric
  if(errorType=='L2'){
    returnValue <- sqrt(sum(covMatrix^2))
  } else if(errorType=='Tr'){
    returnValue <- sum(diag(covMatrix))
  } else if(errorType=='CE'){
    returnValue <- sum(colSums(abs(CE_std)^2)/M)#sum(colSums(abs(CE_std)^2)/M)
  } else{
    stop('Only L2 and Tr error functions implemented')
  }

  returnValue
}

.split_on_NA <- function(vec) {
  is.sep <- is.na(vec)
  split(vec[!is.sep], cumsum(is.sep)[!is.sep])
}

