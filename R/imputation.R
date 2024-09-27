#' Imputation
#'
#' @param data funts object
#' @param method String to indicate method of imputation.
#'  \itemize{
#'    zero: Fill missing with 0
#'    mean_obs: Fill missing with mean of each observation
#'    median_obs: Fill missing with median of each observation
#'    mean_data: Fill missing with mean at time through all data
#'    median_data: Fill missing with median at time through all data
#'    linear: Linearly interpolate missing data
#'  }
#' @param obs_share_data Boolean in linear interpolation that indicates if
#'  data should be shared across observations. For example, is the end of
#'  observation i related to the start of observation i+1. Default is FALSE,
#'  which suggests independence. If true, the distance between the end and
#'  start of observations is taken to be the mean of the intraobs.
#'
#' @return funts object with interpolated missing data
#' @export
#'
#' @examples
#' temp <- data.frame(c(NA,NA,3:9,NA),
#'                    c(NA,rnorm(2),NA, rnorm(6)),
#'                    rnorm(10),
#'                    c(rnorm(4),rep(NA,3), rnorm(3)),
#'                    rep(NA,10),
#'                    c(rnorm(1), rep(NA,9)),
#'                    c(rnorm(9),NA),
#'                    rnorm(10),
#'                    rnorm(10),
#'                    c(NA,NA,3:9, NA))
#' impute(temp, method='mean_obs')
#' impute(temp, method='mean_data')
#' impute(temp, method='linear')
#' impute(temp, method='linear', obs_share_data=TRUE)
impute <- function(data,
                   method = c('zero','mean_obs','median_obs',
                              'mean_data','median_data',
                              'linear','functional'),
                   obs_share_data = FALSE, ...){
  # TODO:: Read https://onlinelibrary-wiley-com.proxy.lib.uwaterloo.ca/doi/pdf/10.1002%2Fsta4.331
  #   Modern multiple imputation with functional data
  data <- .check_data(data,check.na = FALSE)
  method <- .verify_input(method,
                          c('zero','mean_obs','median_obs', 'mean_data',
                            'median_data', 'linear', 'functional') )

  switch(method,
         zero = {
           imputed_zero = replace(data$data, is.na(data$data), 0)
         },
         mean_obs = {
           apply(data$data,MARGIN = 2,function(x){
             tmp <- replace(x, is.na(x), mean(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           })
         },
         median_obs = {
           apply(data$data,MARGIN = 2,function(x){
             tmp <- replace(x, is.na(x), median(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           })
         },
         mean_data = {
           t(apply(data$data,MARGIN = 1,function(x){
             tmp <- replace(x, is.na(x), mean(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           }))
         },
         median_data = {
           t(apply(data$data,MARGIN = 1,function(x){
             tmp <- replace(x, is.na(x), median(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           }))
         },
         linear={
           .linear_imputatation(data, obs_share_data)
         },
         functional={
           na_row <- rowSums(is.na(data$data))
           data_tmp <- stats::na.omit(data$data)

           data_fd <- fda::eval.fd(evalarg = data$intraobs,
                                   fda::Data2fd(
                                     argvals = data$intraobs[na_row==0], y = data_tmp,
                                     basisobj = fda::create.bspline.basis(
                                       rangeval = range(data$intraobs[na_row==0])) )
           )
           data$data <- ifelse(is.na(data$data),data_fd,data$data)

           data
         },
         # pca={
         #
         # },
         {
           stop(paste0('Incorrect method selected. See documentation.'),
                call. = FALSE)
         })
}


#' Linear Imputation
#'
#' @param data funts object
#' @param obs_share_data Boolean to indicate if the end of an observation
#'  should inform the next. Default is FALSE.
#'
#' @return funts object with imputed values
.linear_imputatation <- function(data, obs_share_data = FALSE) {
  ## TODO:: STUDY,  Modern multiple imputation with functional data
  data <- .check_data(data,check.na = FALSE)
  data_fill <- data$data
  n <- ncol(data$data)
  r <- nrow(data$data)

  for(i in 1:n){
    # Skip if nothing missing
    if(sum(is.na(data_fill[,i]))==0) next

    if(obs_share_data){
      # If i=1, i=n, or else
      if(i==1){
        tmp_y <- as.numeric(c(NA,data_fill[,i],data$data[1,i+1]))
      }else if(i==n){
        tmp_y <- as.numeric(c(data$data[r,i-1],data_fill[,i],NA))
      }else{
        tmp_y <- as.numeric(c(data$data[r,i-1],
                              data_fill[,i],
                              data$data[1,i+1]))
      }

      tmp_x <- as.numeric(c(min( data$intraobs)-mean(diff(data$intraobs)),
                            data$intraobs,
                            max(data$intraobs)+mean(diff(data$intraobs))))
    }else{
      tmp_y <- as.numeric(data_fill[,i])
      tmp_x <- as.numeric(data$intraobs)
    }

    tmp_approx <- tryCatch({
      (stats::approxfun(tmp_x,tmp_y,method = 'linear'))(tmp_x)
    },error=function(e){
      data_fill[,i]
    })

    if(obs_share_data){
      data_fill[,i] <- tmp_approx[-c(1,length(tmp_approx))]
    } else{
      data_fill[,i] <- tmp_approx
    }

    ## Fix remaining NA
    #   Possible from starting/ending with NA
    if(sum(is.na(data_fill[,i]))>0){
      # Skip if not sharing but all NA
      if(sum(is.na(data_fill[,i]))==r && !obs_share_data) next

      # Take Average slope down/up for first/last, or not sharing
      if(i==1 || i==n || !obs_share_data){
        data_i <- data_fill[,i]
        obs_st <- min(which(!is.na(data_i)))
        obs_en <- max(which(!is.na(data_i)))
        slope <-  (data_i[obs_en]-data_i[obs_st]) /
          (data$intraobs[obs_en]-data$intraobs[obs_st])

        # If only 1 observation, share it
        if(is.nan(slope)) {
          data_fill[, i] <- data_fill[obs_st,i]
          next
        }

        # Get missing, sort so that I go inside to outside
        obs_missing <- which(is.na(data_i))
        obs_missing <-
          obs_missing[order(abs((obs_missing-obs_en)/2))]

        for(j in obs_missing){
          # Up or down
          if(j < obs_st){
            data_fill[j, i] <- data_fill[j+1, i]-slope*diff(data$intraobs)[j]
          }else{
            data_fill[j, i] <- data_fill[j-1, i]+slope*diff(data$intraobs)[j-1]
          }
        }
      } else{
        data_i <- data_fill[,i]
        obs_st <- min(which(!is.na(data_i)))
        obs1_st <- max(which(!is.na(data_fill[,i-1])))
        obs_en <- max(which(!is.na(data_i)))
        obs1_en <- min(which(!is.na(data_fill[,i+1])))

        low_x <- (data$intraobs[1] -
                    data$intraobs[ length(data$intraobs)-obs1_st+1 ]) -
          mean(diff(data$intraobs))
        slope_down <- (data_i[obs_st]-data_fill[obs1_st,i-1]) /
          (data$intraobs[obs_st] - low_x)
        up_x <- (data$intraobs[length(data$intraobs)] +
                   data$intraobs[ obs1_en ]) +
          mean(diff(data$intraobs))
        slope_up <- (data_fill[obs1_en,i+1]-data_i[obs_en]) /
          ( up_x - data$intraobs[obs_en] )

        # Get missing, sort so that I go inside to outside
        obs_missing <- which(is.na(data_i))
        obs_missing <-
          obs_missing[order(abs((obs_missing-obs_en)/2))]
          # obs_missing[order(abs(obs_missing-r/2))]

        for(j in obs_missing){
          # Up or down
          if(j < obs_st){
            data_fill[j, i] <- data_fill[j+1, i]-slope_down*diff(data$intraobs)[j]
          }else{
            data_fill[j, i] <- data_fill[j-1, i]+slope_up*diff(data$intraobs)[j-1]
          }
        }

      }
    }
  }

  data$data <- data_fill




  data
}


#' #' Impute data using functional data
#' #'
#' #' This function imputes missing data by fitting a basis.
#' #'
#' #' @param data Numeric data.frame with rows for evaluated values and columns
#' #'    indicating FD
#' #' @param evalPts Numeric vector indicating the evaluated points for each row
#' #'    of data
#' #' @param basis basis object from fda to fit data
#' #' @param ... (Optional) Additional parameters to pass into fda
#' #'
#' #' @return Numeric data.frame with rows for evaluated values and columns
#' #'    indicating FD and no missing values
#' #' @export
#' #'
#' #' @examples
#' #' data_missing <- data.frame(
#' #'   "FD1" = c(1:8, NA, 10) + rnorm(10),
#' #'   "FD2" = 21:30 + rnorm(10),
#' #'   "FD3" = c(1:3, NA, rep(6, 3), 8, rep(NA, 2)) + rnorm(10)
#' #' )
#' #' evalPts <- c(1:10)
#' #' functional_imputation(data_missing, evalPts)
#' functional_imputation <- function(data, evalPts = 1:nrow(data),
#'                                   basis = fda::create.bspline.basis(
#'                                     nbasis = 21,
#'                                     rangeval = c(min(evalPts), max(evalPts))
#'                                   ),
#'                                   ...) {
#'   data_evaled <- matrix(nrow = length(evalPts), ncol = ncol(data))
#'   for (i in 1:ncol(data)) {
#'     data_tmp <- data.frame(
#'       "evalPts" = evalPts,
#'       "y" = data[, i]
#'     )
#'     data_tmp <- stats::na.omit(data_tmp)
#'     tmp <- fda::Data2fd(
#'       argvals = data_tmp[, 1], y = as.matrix(data_tmp[, 2]),
#'       basisobj = basis, ...
#'     )
#'     data_evaled[, i] <- fda::eval.fd(evalPts, tmp)
#'   }
#'
#'   data_evaled
#' }
#'
#'
#' #' Impute data using linear function
#' #'
#' #' This function imputing missing data using a line.
#' #'
#' #' @param data Numeric data.frame with rows for evaluated values and columns
#' #'    indicating FD
#' #' @param evalPts Numeric vector indicating the evaluated points for each row
#' #'    of data
#' #' @param use.prev.curve (Optional) Boolean indicating if the last functional
#' #'    observation should be used to impute the first value if needed. Default is
#' #'    FALSE.
#' #'
#' #' @return Numeric data.frame with rows for evaluated values and columns
#' #'    indicating FD and no missing values
#' #' @export
#' #'
#' #' @examples
#' #' data_missing <- data.frame(
#' #'   "FD1" = c(1:8, NA, 10) + rnorm(10),
#' #'   "FD2" = 21:30 + rnorm(10),
#' #'   "FD3" = c(1:3, NA, rep(6, 3), 8, rep(NA, 2)) + rnorm(10)
#' #' )
#' #' evalPts <- c(1:10)
#' #' linear_imputatation(data_missing, evalPts)
#' linear_imputatation <- function(data, evalPts = 1:nrow(data),
#'                                 use.prev.curve = FALSE) {
#'   # Look at all FDs
#'   for (i in 1:ncol(data)) {
#'     # If missing value
#'     if (sum(is.na(data[, i])) > 0) {
#'       # Set values if none in column and want to use others
#'       #   Assume filled last and will try to use future, but NA is possible
#'       if (sum(is.na(data[, i])) == nrow(data) &&
#'         use.prev.curve) {
#'         if (i > 1) data[1, i] <- data[nrow(data), i - 1]
#'         if (i < ncol(data)) data[nrow(data), i] <- data[1, i + 1]
#'       }
#'
#'       # Set starting
#'       prevInfo <- c(data[1, i], 1)
#'       nextInfo <- c(data[nrow(data), i], nrow(data))
#'
#'       # Check first
#'       if (is.na(prevInfo[1])) {
#'         for (j in 2:nrow(data)) {
#'           if (!is.na(data[j, i])) {
#'             prevInfo <- c(data[j, i], j)
#'             for (k in 1:j) {
#'               data[k, i] <- prevInfo[1]
#'             }
#'             break
#'           }
#'         }
#'       }
#'       # Check last
#'       if (is.na(nextInfo[1])) {
#'         for (j in (nrow(data) - 1):1) {
#'           if (!is.na(data[j, i])) {
#'             nextInfo <- c(data[j, i], j)
#'             for (k in nrow(data):j) {
#'               data[k, i] <- nextInfo[1]
#'             }
#'             break
#'           }
#'         }
#'       }
#'
#'       # Fix Middle
#'       st <- prevInfo[2]
#'       en <- nextInfo[2]
#'       fill <- FALSE
#'
#'       ## FIX HERE
#'       # if(!is.numeric(st) && !is.numeric(en)){
#'       #   warning(paste0('Error: Column ',j,' has no data in it.',
#'       #                  ' It is entirely dropped from data'))
#'       # } else if()
#'
#'       for (j in (st + 1):en) {
#'         # Check if value needs to be interpolated
#'         if (is.na(data[j, i]) && !fill) {
#'           prevInfo <- c(data[j - 1, i], j - 1)
#'           fill <- TRUE
#'         }
#'         # Interpolate values if possible
#'         if (!is.na(data[j, i]) && fill) {
#'           nextInfo <- c(data[j, i], j)
#'           fill <- FALSE
#'           for (k in (prevInfo[2] + 1):(nextInfo[2] - 1)) {
#'             x1 <- evalPts[prevInfo[2]]
#'             x2 <- evalPts[k]
#'             x3 <- evalPts[nextInfo[2]]
#'             y1 <- prevInfo[1]
#'             y3 <- nextInfo[1]
#'
#'
#'             data[k, i] <- (x2 - x1) * (y3 - y1) / (x3 - x1) + y1
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   data
#' }
