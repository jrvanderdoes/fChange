#' Principal Component Exploration
#'
#' Computes the principal component projection and re-combined data with a myriad
#'  of figures to understand the projection process.
#'
#' @inheritParams pca
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#'
#' @return List with the following elements:
#' \itemize{
#'  \item figures: List with figures on all components, summaries, reconstructed
#'    values, and residuals.
#'  \item reconstruction: Reconstructed value from PCs.
#'  \item residuals: Difference between true and reconstruction.
#' }
#' @export
#'
#' @seealso [pca.dfts()]
#'
#' @examples
#' results <- pca_examination(electricity, TVE=0.9)
pca_examination <- function(X, TVE=0.95){
  X <- dfts(X)

  pc_data <- pca(X,TVE=TVE)
  num_pcs <- length(pc_data$sdev)
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
      paste0('Reconstruction With ', num_pcs,
             'Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  residuals <- X$data - reconstruction
  plot_residuals <-
    .plot_rainbow(dfts(residuals)) +
    ggplot2::ggtitle(
      paste0('Residuals from Reconstruction With ', num_pcs,
             ' Principal Components') ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  list(figures =list(all = plots,
                     summary = plots1,
                     recontruction = plot_reconstruction,
                     residuals = plot_residuals),
       reconstruction = dfts(reconstruction,
                             name = paste0('PCA reconstruction of ',X$name),
                             labels = X$labels,
                             fparam = X$fparam),
       residuals = dfts(residuals,
                        name = paste0('Residuals from PCA reconstruction of ',X$name),
                        labels = X$labels,
                        fparam = X$fparam))
}


#' Functional PCA Components
#'
#' Evaluates the given pca components as a function.
#'
#' @param pca PCA object from `pca()`
#' @param components Numeric for the components of interest. Can be a single
#'  numeric to examine that component or multiple to examine the combined result
#'
#' @returns A dfts object of the pca component(s)
#' @export
#'
#' @examples
#' tmp <- pca(electricity, TVE=0.1)
#' pca_components(tmp, components=1)
pca_components <- function(pca, components=1:length(pca$sdev)){
  pca_coef <- pca$x[,components,drop=FALSE] %*% t(pca$rotation)[components,,drop=FALSE]

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

  dfts(t(val))
}


#' Projection-based functional data model
#'
#' Model and forecast functional data using a Hyndman and Ullah projection-based
#'  model.
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param TVE Numeric in \[0,1\] for the total variance explained to select
#'  number of PCA components to use to model the data.
#' @param forecast.model String to indicate method to model components, either
#'  "ets" or "arima".
#' @param n.ahead Number of observations to forecast.
#' @param alpha Significance in \[0,1\] for intervals when forecasting.
#' @param check.cp Boolean which indicates if the errors should be checked for
#'  change points to change forecasts and plots.
#' @param sim.bounds Boolean if the confidence bounds should be simulated or
#'  computed using the covariance.
#' @param M Numeric for the number of iterations used to simulated confidence
#'  bounds when sim.bounds is TRUE.
#' @param transformation Argument that specifies any transformations. Currently only
#'  NULL (no transformation) and 'log' (logarithmic) are acceptable.
#' @param ... Additional information to pass into pca, change (if
#'  \code{check.cp=TRUE}), and plot.
#'
#' @return List with the following elements:
#' \itemize{
#'  \item data: List with data information.
#'  \item plots: List with various plots.
#'  \item residuals: dfts object for residuals from the fit.
#'  \item changes: vector of any changes when using \code{detect.cp}.
#'  \item parameters: List with fit parameters like pcs, TVE, model, and n.ahead.
#' }
#' @export
#'
#' @references Hyndman, R. J., & Shahid Ullah, M. (2007). Robust forecasting
#'  of mortality and fertility rates: A functional data approach. Computational
#'  Statistics & Data Analysis, 51(10), 4942-4956.
#'  https://doi.org/10.1016/j.csda.2006.07.028
#'
#' @examples
#' result <- projection_model(dfts(electricity$data[,50:100], season=7),
#'  n.ahead=1, TVE=0.1, check.cp=FALSE, sim.bounds=FALSE)
projection_model <- function(X, TVE = 0.95, forecast.model=c('ets','arima'),
                    n.ahead=0, alpha=0.05, check.cp=TRUE,
                    sim.bounds = TRUE, M=1000, transformation=NULL, ...){
  if(!requireNamespace('forecast',quietly = TRUE))
    stop("Please install the 'forecast' package",call. = FALSE)
  # Verify input
  forecast.model <- .verify_input(forecast.model, c('ets','arima'))
  if(TVE >1 || TVE <0){
    TVE <- max(min(TVE,1),0)
    warning('TVE should be in [0,1]. It was automatically rounded to ensure this.')
  }

  if(is.null(transformation)){
    X <- dfts(X)
    transformation <- 'none'
  }else if(transformation=='log'){
    Xmin <- min(X$data)
    X <- log(dfts(X)+Xmin+1)
  } else{
    X <- dfts(X)
    transformation <- 'none'
    warning('Check transformation argument. No transformation used')
  }

  # Project data and model each component
  pc_data <- pca(X, TVE = TVE, ...)

  comps <- comps_fits <- comps_true <- list()
  for (i in 1:ncol(pc_data$x)) {
    if(forecast.model=='ets'){
      comps[[i]] <- forecast::ets(stats::ts(pc_data$x[, i],frequency =X$season))
    }else if(forecast.model=='arima'){
      comps[[i]] <- forecast::auto.arima(stats::ts(pc_data$x[, i],frequency = X$season))
    }

    # Forecast as request
    for_vals <- NULL
    if(n.ahead>0){
      for_vals <- forecast::forecast(comps[[i]], h=n.ahead)
    }
    comps_fits[[i]] <- c(comps[[i]]$fitted, for_vals$mean)
    comps_true[[i]] <- pc_data$x[, i]
  }
  data_fits <- do.call(cbind, comps_fits)
  data_reals <- do.call(cbind, comps_true)

  # Create functional data and its error
  fit <- .pca_reconstruct(pc_data, data_fits)
  errors <- X$data - fit[,1:ncol(X)]

  if(n.ahead>0){
    fit_prep <- dfts(fit, name = 'Fit', fparam=X$fparam)
    data_prep <- dfts(cbind(X$data,fit[,(ncol(X)+1):ncol(fit)]),
                       name = 'Forecast', fparam=X$fparam)
  }else{
    fit_prep <- dfts(fit, labels=X$labels, name='Fit', fparam=X$fparam,)
    data_prep <- dfts(X$data, labels=X$labels,
                       name = 'Forecast', fparam=X$fparam)
  }

  # Check changes
  if(check.cp){
    res <- fchange(X = dfts(errors), type='segmentation', ...)
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
  lower <- upper <- NULL
  if(n.ahead>0){

    lower <- upper <- matrix(nrow=nrow(X),ncol=n.ahead)

    if(sim.bounds){
      x_b <- array(NA, dim=c(nrow(X), n.ahead,M))
      for(m in 1:M){
        sims <- matrix(nrow=n.ahead,ncol=ncol(pc_data$x))
        for(i in 1:ncol(sims)) {
          sims[,i] <- stats::simulate(comps[[i]], future=TRUE, nsim=n.ahead)
        }
        x_b[,,m] <- .pca_reconstruct(pc_data,sims) +
          errors_use$data[,sample(1:ncol(errors_use),size = n.ahead)]
      }
      lower <- apply(x_b, MARGIN = 2, function(x, alpha){
        apply(x, MARGIN = 1, stats::quantile, probs=c(alpha/2))
      }, alpha=alpha)
      upper <- apply(x_b, MARGIN = 2, function(x, alpha){
        apply(x, MARGIN = 1, stats::quantile, probs=c(1-alpha/2))
      }, alpha=alpha)
    }else{
      covs <- autocovariance(dfts(errors_use),lags = 0)
      for(i in 1:n.ahead){
        lower[,i] <- fit[,ncol(X)+i] - stats::qnorm(1-alpha/2) * sqrt(diag(covs))
        upper[,i] <- fit[,ncol(X)+i] + stats::qnorm(1-alpha/2) * sqrt(diag(covs))
      }
    }

    if(transformation=='log'){
      data_prep <- exp(data_prep)-Xmin-1
      fit_prep <- exp(fit_prep)-Xmin-1
      lower <- exp(lower)-Xmin-1
      upper <- exp(upper)-Xmin-1
    }
    plt_for <- .plot_forecast(data_prep, lower, upper, changes=changes, ...)
    plt_for_fit <- .plot_forecast(fit_prep, lower, upper, changes=changes, ...)
  } else{

    if(transformation=='log'){
      data_prep <- exp(data_prep)-Xmin-1
      fit_prep <- exp(fit_prep)-Xmin-1
      lower <- exp(lower)-Xmin-1
      upper <- exp(upper)-Xmin-1
    }
    plt_for <- plot(data_prep, changes=changes, ...)
    plt_for_fit <- plot(fit_prep, changes=changes, ...)
  }

  # Get Component fits
  components_plots <- list()
  x <- y <- NULL
  for(i in 1:ncol(data_fits)){
    tmp_dat <- data_fits[1:ncol(X),i]
    tmp_dat_for <- data_fits[ncol(X)+1:n.ahead,i]

      if(transformation=='log'){
        tmp_dat <- exp(tmp_dat)-Xmin-1
        tmp_dat_for <- exp(tmp_dat_for)-Xmin-1
      }


      plt <-
        ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x=x, y=y),
                           data=cbind('x'=1:ncol(X),'y'=tmp_dat),
                           col='black') +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                       axis.title = ggplot2::element_text(size=22))
    if(n.ahead>0){
      plt <- plt +
        ggplot2::geom_line(ggplot2::aes(x=x, y=y),
                           data=cbind('x'=ncol(X)+1:n.ahead,
                                      'y'=tmp_dat_for),
                           col='red')

    }

    components_plots[[eval(paste('Component',i))]] <- plt
  }

  list(data = list(
         component_model = fit_prep,
         component_true = data_prep,
         residuals = dfts(errors, labels=X$labels, name='Residuals', fparam=X$fparam)
       ),
       plots = list(
         forecast_plot = plt_for,
         fit_plot = plt_for_fit,
         components = components_plots
       ),
       changes = changes,
       parameters = list(
         pcs = length(pc_data$sdev),
         TVE = TVE,
         forecast.model = forecast.model,
         n.ahead = n.ahead,
         lower = lower,
         upper = upper
       )
  )
}


#' PCA reconstruction
#'
#' Reconstruct data from pca
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
