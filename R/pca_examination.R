#' Forecast using principal component analysis
#'
#' TODO:: Add n to forecast ahead
#'
#' @param X funts object
#' @param TVE Numeric for the total variance explained to select number
#'  of components in PCA
#' @param model String to indicate method to model components, either
#'  "ets" or "arima".
#' @param ... Additional parameters to pass into pca
#'
#' @return List with the following elements:
#' \itemize{
#'  \item fit: funts object for fit
#'  \item errors: funts object for errors from fit
#'  \item parameters: list with fit parameters like pcs, TVE, and model
#' }
#' @export
#'
#' @examples
#' pca_fit(funts(electricity))
#' pca_fit(generate_brownian_motion(100))
pca_fit <- function(X, TVE = 0.95, model=c('ets','arima'), ...){
  if(!requireNamespace('forecast',quietly = TRUE))
    stop("Please install the 'forecast' package",call. = FALSE)

  pc_data <- pca(X, TVE = TVE, ...)

  ## Forecast Each
  comps <- comps_fits <- comps_resids <- list()
  model <- c('ets','arima')[min(pmatch(model,c('ets','arima')))]
  for (i in 1:ncol(pc_data$rotation)) {
    if(model=='ets'){
      comps[[i]] <- forecast::ets(stats::ts(pc_data$rotation[, i]))
    }else if(model=='arima'){
      comps[[i]] <- forecast::auto.arima(stats::ts(pc_data$rotation[, i]))
    }
    comps_resids[[i]] <- stats::residuals(comps[[i]])#$residuals
    comps_fits[[i]] <- comps[[i]]$fitted
  }

  data_errors <- do.call(cbind, comps_resids)
  data_fits <- do.call(cbind, comps_fits)

  # Revert Back to fd
  #   t(data_errors+data_fits) == pc_data$x
  # if(!pc_data$scale){
  #   rebuild <- t(t(pc_data$x %*% t(pc_data$rotation)) +
  #                  pc_data$center)
  # }else{
  #   rebuild <- t(t(pca$x %*% t(pca$rotation)) * pca$scale +
  #                  pca$center)
  # }
  fit <- .pca_reconstruct(pc_data, data_fits)
  # errors <- .pca_reconstruct(pc_data, data_errors)
  errors <- X$data - fit
  # true <- .pca_reconstruct(pc_data)

  # orig_coefs <- pc_data$x %*% rebuild
  # eval_fd_vals <- pc_data$scores %*% orig_coefs

  list(fit = funts(fit),
       errors = funts(errors),
       parameters = list(
         pcs = length(pc_data$sdev),
         TVE = TVE,
         model = model
       ))
}


#' PCA reconstruction
#'
#' @param pca PCA object from pca
#' @param new_rot Forecasted rotation
#'
#' @return Data.frame of forecasts
#' @export
#'
#' @examples
#' pc_data <- pca(funts(electricity),TVE=0.8)
#' fit <- .pca_reconstruct(pc_data)
#' plot_fd(fit)
.pca_reconstruct <- function(pca, new_rot = NULL){
  # if(is.null(pcs)) pcs <- seq(pca$sdev)

  if(is.null(new_rot)){
    pca_coef <- pca$x %*% t(pca$rotation)
  } else{
    pca_coef <- pca$x %*% t(new_rot)
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


#' Principal Component Exploration
#'
#' @param X funts object
#' @param order Number of components to investigate
#' @param ... Unused additional inputs
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
#' @references Hyndman, R. J., & Shahid Ullah, M. (2007). Robust forecasting
#'  of mortality and fertility rates: A functional data approach. Computational
#'  Statistics & Data Analysis, 51(10), 4942–4956.
#'  https://doi.org/10.1016/j.csda.2006.07.028
#'
#' @examples
#' results <- pcaExploration(funts(electricity))
pcaExploration <- function(X, order = 3, ...){
  X <- .check_data(X)

  pc_data <- pca(X)

  results <- array(dim=c(dim(X$data),order))
  plots <- plots1 <- list()
  for(i in 1:order){
    results[,,i] <- .pca_mult(pc_data, i)

    plots[[i]] <- rainbow_plot(funts(results[,,i])) +
      ggplot2::ggtitle(paste0('Principal Component ',i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plots1[[i]] <- distribution_plot(funts(results[,,i])) +
      ggplot2::ggtitle(paste0('Principal Component ',i)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }

  # pc_order <- list(sdev = pc_data$sdev[1:order],
  #                  rotation = pc_data$rotation[,1:order],
  #                  center = pc_data$center,
  #                  scale = pc_data$scale,
  #                  x = pc_data$x[,1:order],
  #                  fullpc = pc_data)
  reconstruction <- .pca_reconstruct(pc_data)

  plot_reconstruction <-
    rainbow_plot(funts(reconstruction)) +
    ggplot2::ggtitle(
      paste0('Reconstruction With ', order,
             'Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  residuals <- X$data - reconstruction
  plot_residuals <-
    rainbow_plot(funts(residuals)) +
    ggplot2::ggtitle(
      paste0('Residuals from reconstruction With ', order,
             ' Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  list(Figures =list(all = plots,
                     summary = plots1,
                     recontruction = plot_reconstruction,
                     residuals = plot_residuals),
       reconstruction = reconstruction,
       residuals = residuals)
}


#' Principal Component Multiplication
#'
#' @param pca PCA object from pca
#' @param pc_num PC of interest
#'
#' @return pca coefficient from the pca and its rotation
#' @export
#'
#' @examples
#' vals <- pca(funts(electricity))
#' res <- .pca_mult(vals,1)
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
