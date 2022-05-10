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
#' results <- elbow_method(data_KL,
#'                   test_statistic_function=compute_Mn, fn=compute_Mn,
#'                   cutoff_function = generalized_resampling,
#'                   trim_function = trim_function,
#'                   alpha = 0.05, M=1000)
#' print(results[[2]])
#' plot_fd(data_KL, CPs=results[[1]]$CP[1:2])
elbow_method <- function(data, test_statistic_function,
                         cutoff_function, trim_function,
                         alpha=0.05, errorType = 'L2', ... ){
  # Setup
  n <- ncol(data)
  return_data <- data.frame('CP'=NA,
                            'Var'=NA)

  # Trim & stopping criteria
  trim_amt <- trim_function(data, ...)
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
        test_stat[k] <- test_statistic_function(data[, prev2currCP],
                                                k-min(prev2currCP), ...)
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
        test_stat[k] <- test_statistic_function(data[, curr2nextCP],
                                                k-min(curr2nextCP), ...)
      }
    }

    if(sum(is.na(test_stat))==n)
      break

    tmp <- c(which.max(test_stat),
             return_data$CP)
    tmp <- tmp[order(tmp)]
    return_data[i,] <- c(which.max(test_stat),
                         .compute_total_var(data, tmp, errorType))
  }

  var_plot <-
    ggplot(return_data) +
    geom_point(aes(x=1:length(CP), y=Var)) +
    geom_line(aes(x=1:length(CP), y=Var)) +
    theme_bw()

  return(list(return_data, var_plot))
}

.compute_total_var <- function(data, CPs, errorType='L2'){
  # Setup
  data_std <- data
  CPs <- c(0, CPs, ncol(data))

  # Loop
  for(i in 2:length(CPs)){
    indices <- (CPs[i-1]+1):CPs[i]
    CP_mean <- rowMeans(data[,indices])

    for(k in indices){
      data_std[,k] <- data_std[,k]-CP_mean
    }
  }

  ## Cannot use cov because it removes the mean from individual FDs, thus
  #     making each change point near equally effective
  covMatrix <- matrix(NA,ncol=ncol(data), nrow = ncol(data))
  for(i in 1:ncol(data)){
    for(j in i:ncol(data)){
      tmp <- rep(NA,nrow(data))
      for(k in 1:nrow(data)){
        tmp[k] <- data_std[k,i] * data_std[k,j]
      }
      tmp <-1/(nrow(data)-1)*tmp
      covMatrix[i,j] <- covMatrix[j,i] <- sum(tmp)
    }
  }

  returnValue <- 0
  if(errorType=='L2'){
    returnValue <- sqrt(sum(covMatrix^2))
  } else if(errorType=='Tr'){
    returnValue <- sum(diag(covMatrix))
  } else{
    stop('Only L2 and Tr error functions implemented')
  }

  returnValue
}
