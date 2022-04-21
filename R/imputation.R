
#' Impute data using functional data
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param evalPts Numeric vector indicating the evaluated points for each row
#'    of data
#' @param basis basis object from fda to fit data
#' @param ... (Optional) Additional parameters to pass into fda
#'
#' @return Numeric data.frame with rows for evaluated values and columns
#'    indicating FD and no missing values
#' @export
#'
#' @examples
#' data_missing <- data.frame('FD1'=c(1:8,NA,10)+rnorm(10),
#'                            'FD2'=21:30+rnorm(10),
#'                            'FD3'=c(1:3,NA,rep(6,3),8,rep(NA,2))+rnorm(10))
#' evalPts <- c(1:10)
#' functional_imputation(data_missing, evalPts)
functional_imputation <- function(data, evalPts = 1:nrow(data),
                            basis=fda::create.bspline.basis(nbasis = 21,
                                      rangeval = c(min(evalPts),max(evalPts))),
                            ...){

  data_evaled <- matrix(nrow=length(evalPts), ncol=ncol(data))
  for(i in 1:ncol(data)){
    data_tmp <- data.frame('evalPts'=evalPts,
                           'y'=data[,i])
    data_tmp <- na.omit(data_tmp)
    tmp <- fda::Data2fd(argvals = data_tmp[,1], y=as.matrix(data_tmp[,2]),
                        basisobj = basis,...)
    data_evaled[,i] <- fda::eval.fd(evalPts,tmp)
  }

  data_evaled
}


#' Impute data using linear function
#'
#' @param data Numeric data.frame with rows for evaluated values and columns
#'    indicating FD
#' @param evalPts Numeric vector indicating the evaluated points for each row
#'    of data
#'
#' @return Numeric data.frame with rows for evaluated values and columns
#'    indicating FD and no missing values
#' @export
#'
#' @examples
#' data_missing <- data.frame('FD1'=c(1:8,NA,10)+rnorm(10),
#'                            'FD2'=21:30+rnorm(10),
#'                            'FD3'=c(1:3,NA,rep(6,3),8,rep(NA,2))+rnorm(10))
#' evalPts <- c(1:10)
#' linear_imputatation(data_missing, evalPts)
linear_imputatation <- function(data, evalPts=1:nrow(data)){

  # Look at all FDs
  for(i in 1:ncol(data)){
    # If missing value
    if(sum(is.na(data[,i]))>0){
      # Set starting
      prevInfo <- c(data[1,i], 1)
      nextInfo <- c(data[nrow(data),i],nrow(data))

      # Check first
      if(is.na(prevInfo[1])){
        for(j in 2:nrow(data)){
          if(!is.na(data[j,i])){
            prevInfo <- c(data[j,i], j)
            for(k in 1:j){
              data[k,i] <- prevInfo[1]
            }
            break
          }
        }
      }
      # Check last
      if(is.na(nextInfo[1])){
        for(j in (nrow(data)-1):1){
          if(!is.na(data[j,i])){
            nextInfo <- c(data[j,i], j)
            for(k in nrow(data):j){
              data[k,i] <- nextInfo[1]
            }
            break
          }
        }
      }

      # Fix Middle
      st <- prevInfo[2]
      en <- nextInfo[2]
      fill <- F

      for(j in (st+1):en){
        # Check if value needs to be interpolated
        if(is.na(data[j,i]) && !fill){
          prevInfo <- c(data[j-1,i], j-1)
          fill <- T
        }
        # Interpolate values if possible
        if(!is.na(data[j,i]) && fill){
          nextInfo <- c(data[j,i], j)
          fill <- F
          for(k in (prevInfo[2]+1):(nextInfo[2]-1)){
            x1 <- evalPts[prevInfo[2]]
            x2 <- evalPts[k]
            x3 <- evalPts[nextInfo[2]]
            y1 <- prevInfo[1]
            y3 <- nextInfo[1]


            data[k, i] <- (x2-x1)*(y3-y1)/(x3-x1) + y1
          }
        }
      }

    }
  }

  data
}
