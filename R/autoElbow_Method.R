
#' Auto Elbow Method
#'
#' Method to determine the number of change points using the elbow method.
#'  Change points are selected based on the max \eqn{M_n} statistic (see
#'  `compute_Mn()`) in sub-segments, selecting the one that maximized total
#'  variance explained.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param trim_function_method Function taking data as an argument and returning a
#'     numeric value indicating how much should be trimmed on each end
#' @param newVarExplain XXXXXXXXX
#' @param maxCP XXXXXXXXX
#' @param errorType String of 'L2' or 'Tr' indicating the error function to use
#' @param demean XXXXXXXXX
#' @param ... Additional parameters to pass into the `trim_function_method()` or
#' `compute_Mn()`.
#'
#' @return list with element 1 the data frame with the change point location and
#'     element a ggplot of variance as a function of CPs
#' @export
#'
#' @examples
#' results <- .autoElbow_method(data=electricity)
.autoElbow_method <- function(data,
                              trim_function_method=trim_function,
                              newVarExplain = 0.05, maxCP = 10,
                              errorType = 'L2', demean=TRUE, ... ){
  # Params that could be entered, but only one option for now
  test_stat_method <- compute_Mn
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


  # Run First
  cat(paste0('Find CP (iters): 1'))

  test_stats <- test_stat_method(data, ...)$allValues
  loc <- which.max(test_stats[nStart:nEnd])+nStart-1

  return_data[1,] <- c(loc, .compute_total_var(data,loc, errorType), 1)

  # Search until max CP or below varChange
  #   Non-recursive as no absolute cutoff and only want to certain number
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
    trim_amt <- trim_function_method(data[,prev2currCP])#, ...)
    test_stats[prev2currCP] <- test_stat_method(data[,prev2currCP])$allValues#, ...)$allValues
    # Make NA if range too small (trimmed)
    test_stats[c(min(prev2currCP) + 1:trim_amt - 1,
                 max(prev2currCP) - 1:trim_amt + 1)] <- NA

    # Next
    curr2nextCP <- (CPs[prev_CP_loc]+1):CPs[prev_CP_loc+1]
    trim_amt <- trim_function_method(data[,curr2nextCP])#, ...)
    test_stats[curr2nextCP] <- test_stat_method(data[,curr2nextCP])$allValues#, ...)$allValues
    # Make NA if range too small (trimmed)
    test_stats[c(min(curr2nextCP) + 1:trim_amt - 1,
                 max(curr2nextCP) - 1:trim_amt + 1)] <- NA

    # If there is no possibilities after trimming
    if(sum(is.na(test_stats))==n) break

    ## Get total variance for each potential CP on the segments
    #   Can do less repeating here eventually
    data_segments <- .split_on_NA(test_stats)

    return_data_tmp <- data.frame('CP'=NA,'Var'=NA)
    for(k in 1:length(data_segments)){
      # Find max test stat on interval
      value_max <- max(data_segments[[k]])

      # Get CP and total variance with full data
      section_max <- which(test_stats==value_max)
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

#' Demean Data
#'
#' This (internal) function de-means the data based on some given change points.
#'
#' @inheritParams .autoElbow_method
#' @param CPs Vector of numerics (integers) indicating the location of the
#'  changes in the data.
#'
#' @return Data.frame of the form given in data, but with the means removed.
#'
#' @noRd
.demean <- function(data, CPs){
  CPs <- c(0,CPs[order(CPs)], ncol(data))
  data_demean <- data

  for(i in 1:(length(CPs)-1)){
    data_demean[,(CPs[i]+1):CPs[i+1]] <-
      data_demean[,(CPs[i]+1):CPs[i+1]] - rowMeans(data_demean[,(CPs[i]+1):CPs[i+1]])
  }

  data_demean
}
