
#' Auto Elbow Method
#'
#' Method to determine the number of change points using the elbow method. Note,
#'     cascading change points are not considered to allow for every possible
#'     number of change points to be selectable.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param test_stat_method XXXXX
#' @param trim_function_method Function taking data as an argument and returning a
#'     numeric value indicating how much should be trimmed on each end
#' @param newVarExplain XXXXXXXXX
#' @param maxCP XXXXXXXXX
#' @param errorType String of 'L2' or 'Tr' indicating the error function to use
#' @param demean XXXXXXXXX
#' @param ... Additional parameters to pass into the respective functions
#'
#' @return list with element 1 the data frame with the change point location and
#'     element a ggplot of variance as a function of CPs
#' @export
#'
#' @examples
#' \dontrun{
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
#' results <- .autoElbow_method(data=data_KL,
#'                               test_stat_method=compute_Tn)
#' }
.autoElbow_method <- function(data,
                              test_stat_method=compute_Mn,
                              trim_function_method=trim_function,
                              newVarExplain = 0.05, maxCP = 10,
                              errorType = 'L2', demean=TRUE, ... ){
  # Setup
  n <- ncol(data)
  return_data <- data.frame('CP'=NA,
                            'Var'=NA,
                            'Change'=NA)

  # Trim & stopping criteria
  trim_amt <- trim_function_method(data, ...)
  nStart <- 1+trim_amt
  nEnd <- ncol(as.data.frame(data))-trim_amt
  if(nStart> nEnd) return()

  cat(paste0('Find CP (iters): 1'))

  # Run First
  test_stat <- rep(NA, n)
  for(k in nStart:nEnd){
    test_stat[k] <- test_stat_method(data, k, ...)
  }
  return_data[1,] <- c(which.max(test_stat),
                       .compute_total_var(data, which.max(test_stat),
                                          errorType),1)

  # Search until max CP or below varChange
  i <- 1
  # Change second or redefine base change value
  while(nrow(return_data)<=min(n-1,maxCP) &&
        return_data[i,'Change']>newVarExplain){
    # Setup
    i <- i + 1
    cat(paste0(', ',i))
    CPs <- c(0,return_data$CP, n)
    CPs <- CPs[order(CPs)]
    prev_CP <- return_data[nrow(return_data),'CP']
    prev_CP_loc <- which(CPs==prev_CP)

    # Explanation:
    #   Rather than searching all intervals, only need to check new
    #   This will be to left and right of previous split
    #   So we search each segment, after appropriately trimming

    # Previous
    prev2currCP <- (CPs[prev_CP_loc-1]+1):(CPs[prev_CP_loc])
    trim_amt <- trim_function_method(data[,prev2currCP], ...)
    for(k in prev2currCP){
      if(k < min(prev2currCP)+trim_amt ||
         k > max(prev2currCP)-trim_amt){
        test_stat[k] <- NA
      }
      else{
        test_stat[k] <- test_stat_method(
          as.data.frame(data[, prev2currCP]),
          k-min(prev2currCP)+1, ...)
      }
    }
    # Next
    curr2nextCP <- (CPs[prev_CP_loc]+1):CPs[prev_CP_loc+1]
    trim_amt <- trim_function_method(data[,curr2nextCP], ...)
    for(k in curr2nextCP){
      if(k < min(curr2nextCP)+trim_amt ||
         k > max(curr2nextCP)-trim_amt){
        test_stat[k] <- NA
      }
      else{
        test_stat[k] <- test_stat_method(
          as.data.frame(data[, curr2nextCP]),
          k-min(curr2nextCP)+1, ...)
      }
    }

    # If there is no possibilities after trimming
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
    return_data[i,] <- c(return_data_tmp[which.min(return_data_tmp$Var),],NA)

    return_data[i,'Change'] <-
      abs(1-return_data[i,'Var']/return_data[i-1,'Var'])
  }
  cat(paste0('\n'))

  return(list('CPs'=return_data$CP,
              'data'=data,
              'data_demean'=.demean(data,return_data$CP)))
}

#' Title
#'
#' @param data XXXXX
#' @param CPs XXXXX
#'
#' @return XXXXX
#'
#' @examples
#' # XXXXX
.demean <- function(data, CPs){
  CPs <- c(0,CPs[order(CPs)], ncol(data))
  data_demean <- data

  for(i in 1:(length(CPs)-1)){
    data_demean[,(CPs[i]+1):CPs[i+1]] <-
      data_demean[,(CPs[i]+1):CPs[i+1]] - rowMeans(data_demean[,(CPs[i]+1):CPs[i+1]])
  }

  data_demean
}
