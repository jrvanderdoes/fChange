
#' Complete Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. Change points are recursively found until no
#'     more change points are detected.
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
#' @param final_verify (Optional) Boolean value
#'
#' Indicates if a final pass looking at sequences with only one change point
#'     should be conducted to verify results. Note, this may modify existing
#'     locations of change points, potentially to less accurate locations.
#' @param silent (Optional) Boolean value
#'
#' Indicates if useful output should be silenced. Default F shows output.
#' @param ... Additional arguments passed into test_statistic_function
#'
#' @return A list of numeric values indicating change points  (if exists),
#'     NA otherwise
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(200),
#'     eigsList = list(c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0),
#'     distsArray = c('Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0))
#' complete_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,1),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' complete_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#'
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5),
#'                     c(3,2)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=2, norder=2)),
#'     meansList = c(-1,1,1),
#'     distsArray = c('Normal','Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0,0))
#' complete_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
complete_binary_segmentation <- function(data,
                                         test_statistic_function,
                                         cutoff_function,
                                         trim_function,
                                         alpha = 0.05,
                                         final_verify = T,
                                         silent = F,
                                         ... ){

  # Get change points -- Will not include first or last
  CPsVals <- .detectChangePoints(data=data,
                                test_statistic_function=test_statistic_function,
                                cutoff_function=cutoff_function,
                                trim_function=trim_function,
                                alpha=alpha,
                                silent = silent,
                                ... )

  # Verify as desired
  if(final_verify){

    CPsVals <- .verify(CPsVals=CPsVals, data=data,
                      test_statistic_function=test_statistic_function,
                      cutoff_function=cutoff_function,
                      trim_function=trim_function,
                      alpha=alpha,
                      silent=silent,
                      ...)
  }

  return(CPsVals)
}


#' Single Binary Segmentation
#'
#' This function implements traditional binary segmentation on functional data
#'     for general functions. At most one change point is detected.
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param test_statistic_function Function with the first argument being data
#'     and the second argument optional argument for candidate change points.
#'     Additional arguments passed in via ... . Return a single numeric value.
#' @param cutoff_function Function with first argument being data and the second
#'     argument being alpha. No other arguments given. Return single numeric
#'     value.
#' @param trim_function Function taking data as an argument and returning a
#'     numeric value indicating how much should be trimmed on each end
#' @param alpha Numeric value in [0,1] indicating the significance for
#'     cutoff_function.
#' @param ... Additional arguments passed into test_statistic_function
#'
#' @return A numeric value indicating the cutoff location (if exists),
#'     NA otherwise
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(200),
#'     eigsList = list(c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(0),
#'     distsArray = c('Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0))
#' single_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#' single_binary_segmentation(data_KL, compute_Mn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#'
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(100,100),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,1),
#'     distsArray = c('Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0))
#' single_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#' single_binary_segmentation(data_KL, compute_Mn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
single_binary_segmentation <- function(data, test_statistic_function,
                                       cutoff_function,
                                       trim_function,
                                       alpha=0.05, include_value=F,
                                       ... ){
  # Trim & stopping criteria
  trim_amt <- trim_function(data, ...)
  nStart <- 1+trim_amt
  nEnd <- ncol(as.data.frame(data))-trim_amt
  if(nStart> nEnd) ifelse(include_value,return(c(NA,NA)),return(NA))

  # Find test statistic at every candidate change point
  test_stat <- rep(NA, ncol(data))
  for(k in nStart:nEnd){
    test_stat[k] <- test_statistic_function(data, k, ...)
  }
  test_stat_full <- max(test_stat, na.rm = T)

  # Return index of max change point if larger than cutoff
  return_value <- ifelse(test_stat_full >= cutoff_function(data, alpha, ...),
                which.max(test_stat),
                NA)

  # Add in value
  if(include_value){
    if(!is.na(return_value)){
      return_value <- c(return_value, max(test_stat, na.rm = T))
    } else{
      return_value <- c(return_value, NA)
    }
  }

  return_value
}

#' Wild Binary Segmentation
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param M Numeric value
#'
#' Indicates the number of intervals to examine
#'
#' @param add_full Boolean value
#'
#' Indicates if the entire interval should be added
#'
#' @param block_size (Optional) Numeric value
#'
#' Indicates the minimum block size
#'
#' @param ... Additional parameters to pass into
#'     \code{single_binary_segmentation}
#'
#' @return Numeric values indicating the change points detected
#' @export
#'
#' @examples
#' # Setup Data
#' data_KL <- generate_data_fd(ns = c(12,12,12),
#'     eigsList = list(c(3,2,1,0.5),
#'                     c(3,2,1,0.5),
#'                     c(3,2,1,0.5)),
#'     basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4),
#'                      fda::create.bspline.basis(nbasis=4, norder=4)),
#'     meansList = c(-1,0,1),
#'     distsArray = c('Normal','Normal','Normal'),
#'     evals = seq(0,1,0.05),
#'     kappasArray = c(0,0,0))
#'
#' complete_binary_segmentation(data_KL, compute_Tn, welsh_approximation,
#'     function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
#' wild_binary_segmentation(data=data_KL,
#'     test_statistic_function=compute_Tn,
#'     cutoff_function=welsh_approximation,
#'     trim_function=function(data){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)})
wild_binary_segmentation <- function(data, M=5000, add_full=T, block_size=1,
                                     ...){
  # Setup
  n <- ncol(data)
  cps <- c()
  result <- matrix(ncol=2, nrow = M+add_full)

  # Test
  if(n<=1) return()

  # Run
  if(add_full){
    result[1,] <- single_binary_segmentation(data, include_value=T, ...)
  }

  for(i in add_full+1:M){
    min_pt <- sample(1:(n-block_size+1), 1)
    max_pt <- sample((min_pt+block_size-1):n, 1)

    # Must return location and value
    # This needs to return location and values!
    result[i,] <- single_binary_segmentation(data[,min_pt:max_pt], include_value=T, ...)
  }

  # Select best and continue if reasonable
  if(nrow(na.omit(result))){
    cp_loc <- result[which.max(result[,2]),1]

    cps <- c(wild_binary_segmentation(data[,min_pt:cp_loc],
                 M=M, add_full=add_full,block_size=block_size,
                 ...),
             cp_loc,
             wild_binary_segmentation(data[,(cp_loc+1):max_pt],
                 M=M, add_full=add_full,block_size=block_size,
                 ...) + cp_loc
          )
  }

  cps
}

## This is multiple single_binary_segmentation for complete_binary_segmentation
.detectChangePoints <-  function(data, test_statistic_function,
                                 cutoff_function,
                                 trim_function,
                                 alpha,
                                 addAmt = 0,
                                 silent = F,
                                 ...){

  potential_cp <-
    single_binary_segmentation(data, test_statistic_function,
                               cutoff_function,
                               trim_function,
                               alpha = alpha, ... )

  # No Change Point Detected
  if(is.na(potential_cp)) return()

  # Display progress
  if(!silent)
    cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+ncol(data),' at ',
               addAmt+potential_cp,'): Segment Data and Re-Search\n'))

  return(c(
    .detectChangePoints(data=data[,1:potential_cp],
                       test_statistic_function=test_statistic_function,
                       cutoff_function=cutoff_function,
                       trim_function=trim_function,
                       alpha=alpha,
                       addAmt=addAmt,
                       silent=silent,
                       ...),
    potential_cp + addAmt,
    .detectChangePoints(data=data[,(potential_cp+1):ncol(data)],
                       test_statistic_function=test_statistic_function,
                       cutoff_function=cutoff_function,
                       trim_function=trim_function,
                       alpha=alpha,
                       addAmt=addAmt+potential_cp-1,
                       silent=silent,
                       ...)
  ))
}


## This is verification step in complete_binary_segmentation
.verify <- function(CPsVals, data, test_statistic_function,
                    cutoff_function,
                    trim_function,
                    alpha,
                    silent = F,
                    ...){
  if(!silent) cat('-- Verify Step --\n')

  if(length(CPsVals)>=1){ # If there was a CP detected
    tmp_cps <- c(1,CPsVals, ncol(data))
    newCPVals <- c(1)
    for(i in 2:(length(tmp_cps)-1)){
      tmp <- .detectChangePoints(data = data[,max(newCPVals):tmp_cps[i+1]],
                                test_statistic_function=test_statistic_function,
                                cutoff_function=cutoff_function,
                                trim_function=trim_function,
                                alpha=alpha,
                                silent = silent,
                                addAmt=max(newCPVals)-1,
                                ... )

      newCPVals <- c(newCPVals, tmp)
    }
    CPsVals <- newCPVals[-1]
  } else{
    CPsVals <- .detectChangePoints(data=data,
                                  test_statistic_function=test_statistic_function,
                                  cutoff_function=cutoff_function,
                                  trim_function=trim_function,
                                  alpha=alpha,
                                  silent = silent,
                                  ... )
  }

  # Order
  if(sum(is.na(CPsVals))==length(CPsVals))
    return(NA)
  CPsVals[order(CPsVals)]
}
