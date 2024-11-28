#' Imputation
#'
#' Impute missing values in functional data.
#'
#' @param X Funts object or data easily convertible. See [funts()]
#' @param method String to indicate method of imputation.
#'  \itemize{
#'    \item zero: Fill missing with 0
#'    \item mean_obs: Fill missing with mean of each observation
#'    \item median_obs: Fill missing with median of each observation
#'    \item mean_data: Fill missing with mean at time through all data
#'    \item median_data: Fill missing with median at time through all data
#'    \item linear: Linearly interpolate missing data
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
#'                    c(NA,stats::rnorm(2),NA, stats::rnorm(6)),
#'                    stats::rnorm(10),
#'                    c(stats::rnorm(4),rep(NA,3), stats::rnorm(3)),
#'                    rep(NA,10),
#'                    c(stats::rnorm(1), rep(NA,9)),
#'                    c(stats::rnorm(9),NA),
#'                    stats::rnorm(10),
#'                    stats::rnorm(10),
#'                    c(NA,NA,3:9, NA))
#' impute(temp, method='mean_obs')
#' impute(temp, method='linear', obs_share_data=TRUE)
impute <- function(X,
                   method = c('zero','mean_obs','median_obs',
                              'mean_data','median_data',
                              'linear','functional'),
                   obs_share_data = FALSE){
  # TODO:: Read https://onlinelibrary-wiley-com.proxy.lib.uwaterloo.ca/doi/pdf/10.1002%2Fsta4.331
  #   Modern multiple imputation with functional data
  X <- funts(X, inc.warnings = F)
  method <- .verify_input(method,
                          c('zero','mean_obs','median_obs', 'mean_data',
                            'median_data', 'linear', 'functional') )

  X_imp <- switch(method,
         zero = {
           imputed_zero = replace(X$data, is.na(X$data), 0)
         },
         mean_obs = {
           apply(X$data,MARGIN = 2,function(x){
             tmp <- replace(x, is.na(x), mean(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           })
         },
         median_obs = {
           apply(X$data,MARGIN = 2,function(x){
             tmp <- replace(x, is.na(x), stats::median(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           })
         },
         mean_data = {
           t(apply(X$data,MARGIN = 1,function(x){
             tmp <- replace(x, is.na(x), mean(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           }))
         },
         median_data = {
           t(apply(X$data,MARGIN = 1,function(x){
             tmp <- replace(x, is.na(x), stats::median(x, na.rm = TRUE))
             replace(tmp,is.nan(tmp),NA)
           }))
         },
         linear={
           .linear_imputatation(X, obs_share_data)
         },
         functional={
           if(!requireNamespace('fda',quietly = TRUE))
             stop('Functional not possible without fda package. Install and re-run code', call. = FALSE)

           na_row <- rowSums(is.na(X$data))
           data_tmp <- stats::na.omit(X$data)

           if(length(X$intraobs[na_row==0])==0){
             stop("Need some evaluation points to be non-NA for all observations",call. = FALSE)
           }

           if(na_row[1]>0 || na_row[length(na_row)]>0){
             warning('fda may fail with NAs in first or last evalation points',call. = FALSE)
           }


           data_fd <- fda::eval.fd(evalarg = X$intraobs,
                                   fda::Data2fd(
                                     argvals = X$intraobs[na_row==0], y = data_tmp,
                                     basisobj = fda::create.bspline.basis(
                                       rangeval = range(X$intraobs[na_row==0])) )
           )
           X$data <- ifelse(is.na(X$data),data_fd,X$data)

           X
         },
         # pca={
         #
         # },
         {
           stop(paste0('Incorrect method selected. See documentation.'),
                call. = FALSE)
         })

  funts(X_imp,name=X$name,labels=X$labels, intraobs=X$intraobs, inc.warnings = F)
}


#' Linear Imputation
#'
#' @param data funts object
#' @param obs_share_data Boolean to indicate if the end of an observation
#'  should inform the next. Default is FALSE.
#'
#' @return funts object with imputed values
#'
#' @keywords internal
#' @noRd
.linear_imputatation <- function(data, obs_share_data = FALSE) {
  ## TODO:: STUDY,  Modern multiple imputation with functional data
  data <- funts(data, inc.warnings = FALSE)
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
