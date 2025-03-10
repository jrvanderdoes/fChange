#' Principal Component Exploration
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param TVE Numeric in \[0,1\]. Total variance explained to determine the
#'  number of components to examine
#' @param d.max Numeric. Max number of components to investigate (if TVE
#'  suggests too many).
#'
#' @return List with the following elements:
#' \itemize{
#'  \item Figures: List with figures on all components, summaries, reconstructed
#'    values, and residuals
#'  \item reconstruction: Reconstructed value from pcs
#'  \item residuals: Difference between true and reconstruction
#' }
#' @export
#'
#' @examples
#' results <- pca_examination(electricity)
pca_examination <- function(X, TVE=0.95, d.max = 3){
  X <- dfts(X)

  pc_data <- pca(X,TVE=TVE)
  num_pcs <- min(length(pc_data$sdev),d.max)
  pc_data$sdev <- pc_data$sdev[1:num_pcs,drop=FALSE]
  pc_data$rotation <- pc_data$rotation[,1:num_pcs,drop=FALSE]
  pc_data$x <- pc_data$x[,1:num_pcs,drop=FALSE]


  results <- array(dim=c(dim(X$data),num_pcs))
  plots <- plots1 <- list()
  for(i in 1:num_pcs){
    results[,,i] <- .pca_mult(pc_data, i)

    plots[[i]] <- .plot_rainbow(dfts(results[,,i])) +
      ggplot2::ggtitle(paste0('Principal Component ',i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plots1[[i]] <- .plot_banded(dfts(results[,,i])) +
      ggplot2::ggtitle(paste0('Principal Component ',i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  # pc_order <- list(sdev = pc_data$sdev[1:d.max],
  #                  rotation = pc_data$rotation[,1:d.max],
  #                  center = pc_data$center,
  #                  scale = pc_data$scale,
  #                  x = pc_data$x[,1:d.max],
  #                  fullpc = pc_data)
  reconstruction <- .pca_reconstruct(pc_data)

  plot_reconstruction <-
    .plot_rainbow(dfts(reconstruction)) +
    ggplot2::ggtitle(
      paste0('Reconstruction With ', d.max,
             'Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  residuals <- X$data - reconstruction
  plot_residuals <-
    .plot_rainbow(dfts(residuals)) +
    ggplot2::ggtitle(
      paste0('Residuals from reconstruction With ', d.max,
             ' Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  list(Figures =list(all = plots,
                     summary = plots1,
                     recontruction = plot_reconstruction,
                     residuals = plot_residuals),
       reconstruction = reconstruction,
       residuals = residuals)
}


#' Forecast using principal component analysis
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param TVE Numeric for the total variance explained to select number
#'  of components in PCA
#' @param model String to indicate method to model components, either
#'  "ets" or "arima".
#' @param n.ahead Number of observations to forecast
#' @param alpha Significance in \[0,1\] for intervals when forecasting
#' @param check.cp Boolean which indicates if the errors should be checked for
#'  change points to change forecasts and plots
#' @param frequency Numeric for seasonal frequency when component is made a ts
#'  object
#' @param ... Additional information to pass into pca, change (if check.cp=TRUE),
#'  and plot
#'
#' @return List with the following elements:
#' \itemize{
#'  \item fit: dfts object for fit
#'  \item errors: dfts object for errors from fit
#'  \item parameters: list with fit parameters like pcs, TVE, and model
#' }
#' @export
#'
#' @references Hyndman, R. J., & Shahid Ullah, M. (2007). Robust forecasting
#'  of mortality and fertility rates: A functional data approach. Computational
#'  Statistics & Data Analysis, 51(10), 4942-4956.
#'  https://doi.org/10.1016/j.csda.2006.07.028
#'
#' @examples
#' result <- projection_model(dfts(electricity$data[,50:100]),
#'  n.ahead=1, TVE=0.1, check.cp=FALSE)
projection_model <- function(X, TVE = 0.95, model=c('ets','arima'),
                    n.ahead=0, alpha=0.05, check.cp=TRUE, frequency=1,
                    ...){
  if(!requireNamespace('forecast',quietly = TRUE))
    stop("Please install the 'forecast' package",call. = FALSE)
  # Verify input
  model <- .verify_input(model, c('ets','arima'))
  if(TVE >1 || TVE <0){
    TVE <- max(min(TVE,1),0)
    warning('TVE should be in [0,1]. It was automatically rounded to ensure this.')
  }
  X <- dfts(X)

  # Project data and model each component
  pc_data <- pca(X, TVE = TVE, ...)

  comps_fits <- comps_true <- list()
  for (i in 1:ncol(pc_data$x)) {
    if(model=='ets'){
      comps <- forecast::ets(stats::ts(pc_data$x[, i],frequency =frequency))
    }else if(model=='arima'){
      comps <- forecast::auto.arima(stats::ts(pc_data$x[, i],frequency = frequency))
    }

    # Forecast as request
    for_vals <- NULL
    if(n.ahead>0){
      for_vals <- forecast::forecast(comps, h=n.ahead)
    }
    comps_fits[[i]] <- c(comps$fitted, for_vals$mean)
    comps_true[[i]] <- pc_data$x[, i]
  }
  data_fits <- do.call(cbind, comps_fits)
  data_reals <- do.call(cbind, comps_true)

  # Create functional data and its error
  fit <- .pca_reconstruct(pc_data, data_fits)
  errors <- X$data - fit[,1:ncol(X)]

  if(n.ahead>0){
    fit_prep <- dfts(fit, name = 'Fit', intratime=X$intratime)
    data_prep <- dfts(cbind(X$data,fit[,(ncol(X)+1):ncol(fit)]),
                       name = 'Forecast', intratime=X$intratime)
  }else{
    fit_prep <- dfts(fit, labels=X$labels, name='Fit', intratime=X$intratime,)
    data_prep <- dfts(X$data, labels=X$labels,
                       name = 'Forecast', intratime=X$intratime)
  }

  # Check changes
  if(check.cp){
    res <- change(X = dfts(errors), type='segmentation', ...)
    if(!is.null(res)){
      changes <- res$location
      errors_use <- dfts(errors[,(max(res$location)+1):ncol(X),drop=FALSE])
    }else{
      errors_use <- dfts(errors)
      changes <- NULL
    }
  } else{
    errors_use <- dfts(errors)
    changes <- NULL
  }

  # Confidence intervals for forecast and plots
  if(n.ahead>0){

    covs <- autocovariance(dfts(errors_use),lags = 0)
    lower <- upper <- matrix(nrow=nrow(X),ncol=n.ahead)
    for(i in 1:n.ahead){
      lower[,i] <- fit[,ncol(X)+i] - stats::qnorm(1-alpha/2) * sqrt(diag(covs))
      upper[,i] <- fit[,ncol(X)+i] + stats::qnorm(1-alpha/2) * sqrt(diag(covs))
    }

    plt_for <- .plot_forecast(data_prep, lower, upper, changes=changes, ...)
    plt_for_fit <- .plot_forecast(fit_prep, lower, upper, changes=changes, ...)
  } else{
    plt_for <- plot(data_prep, changes=changes, ...)
    plt_for_fit <- plot(fit_prep, changes=changes, ...)
  }

  list(fit = fit_prep,
       forecast_plot = plt_for,
       fit_plot = plt_for_fit,
       errors = dfts(errors, labels=X$labels, name='Errors', intratime=X$intratime),
       changes = changes,
       component_models = data_fits,
       component_true = data_reals,
       parameters = list(
         pcs = length(pc_data$sdev),
         TVE = TVE,
         model = model,
         n.ahead = n.ahead
       ) )
}


#' PCA reconstruction
#'
#' Recontruct data from pca
#'
#' @param pca PCA object from pca
#' @param new_pca New pca object
#'
#' @return Data.frame of forecasts
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' #pc_data <- pca(electricity,TVE=0.8)
#' #fit <- .pca_reconstruct(pc_data)
.pca_reconstruct <- function(pca, new_pca = NULL){
  # if(is.null(pcs)) pcs <- seq(pca$sdev)

  if(is.null(new_pca)){
    pca_coef <- pca$x %*% t(pca$rotation)
  } else{
    pca_coef <- new_pca %*% t(pca$rotation)
  }

  if(pca$center[1]){
    center_val <- -pca$center
  }else{
    center_val <- FALSE
  }

  if(pca$scale[1]){
    scale_val <- 1/pca$scale
  }else{
    scale_val <- FALSE
  }

  val <- scale(pca_coef, center=center_val, scale=scale_val )
  attr(val,'scaled:center') <- NULL
  attr(val,'scaled:scale') <- NULL #TODO::CHECK

  t(val)
  # if(pca$scale[1] != FALSE){
  #   pca_coef <- scale(pca_coef, center=FALSE, scale=1/pca$scale)
  # }
  # if(pca$center[1] != FALSE){
  #   pca_coef <- scale(pca_coef, center=-pca$center, scale=FALSE)
  # }
  # pca_coef
}


#' Principal Component Multiplication
#'
#' @param pca PCA object from pca
#' @param pc_num PC of interest
#'
#' @return pca coefficient from the pca and its rotation
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' #vals <- pca(electricity)
#' #res <- .pca_mult(vals,1)
.pca_mult <- function(pca, pc_num = 1){
  pca_coef <- pca$x[,pc_num] %*% t(pca$rotation[,pc_num])

  if(pca$center[1]){
    center_val <- -pca$center
  }else{
    center_val <- FALSE
  }

  if(pca$scale[1]){
    scale_val <- 1/pca$scale
  }else{
    scale_val <- FALSE
  }

  val <- scale(pca_coef, center=center_val, scale=scale_val )
  attr(val,'scaled:center') <- NULL
  attr(val,'scaled:scale') <- NULL #TODO::CHECK

  val
}
