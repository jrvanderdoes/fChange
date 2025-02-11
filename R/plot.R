#' Plot Funts
#'
#' A function to plot funts objects in many ways.
#'
#' @param X Funts object or data easily convertible. See [funts()].
#' @param CPs Vector of change points. Can be NULL.
#' @param type Choice of plotting method. Options include: 'spaghetti', 'highdim',
#'  'rainbow','banded','acf', 'pacf', 'summary', 'qq', 'distribution', 'change',
#'  'interval', and'surface'.
#' @param plot_title Title to include on the return plot
#' @param val_axis_title,res_axis_title,FD_axis_title Title for the axis giving
#'  the values (val), the resolution of the intraobs (res), and the functional
#'  observations (FD)
#' @param eye,aspectratio Angle (eye) and ratio (aspectratio) to view 3d plots
#' @param showticklabels Boolean if the tick marks should be shown
#' @param lag.max Max number of lags to consider for ACF/PACF and summary plots
#' @param d.max Max number of dimensions for qq/distribution and summary plots
#' @param alpha Significance level to be used in various plots. Value in \[0,1\]
#' @param TVE Total variance explained used in qq/distribution plots. Value in \[0,1\]
#' @param distribution String of the distribution to compare against in
#'  distribution plot. The string can be anything such that there is a
#'  rdistribution and ddistribution function available. For example "exp",
#'  "gamma". Additional parameters can be passed using ...
#' @param method Method for computing ACF/PACF. Options include 'Welch', 'MC',
#'  and 'Imhof'.
#' @param legend Boolean if legend should be given in qq/distribution plots
#' @param highlight_changes,int.gradual Boolean if changes should be highlighted
#'  or confidence interval be solid gray or gradual colors in change plot
#' @param intervals Information on confidence intervals of changes for change
#'  plot. See [confidence_interval()].
#' @param ... Details for plotting in acf/pacf, summary, or distribution function
#'
#' @returns Plot of varying types
#' @export
#'
#' @references John Fox, & Sanford Weisberg (2019). An R Companion to Applied
#'  Regression. Sage.
#'
#' @examples
#' plt <- plot(funts(electricity))
#' plt <- plot(var(funts(electricity)), type='surface')
plot.funts <- function(x, CPs=NULL,
                       type=c('spaghetti','highdim', 'rainbow','banded',
                              'acf', 'pacf', 'summary', 'qq', 'distribution',
                              'change','interval', 'surface'),
                       plot_title = x$name, val_axis_title = NULL,
                       res_axis_title = NULL, FD_axis_title = NULL,
                       eye = NULL, aspectratio = NULL,
                       showticklabels = TRUE,
                       lag.max=20, d.max=2,
                       alpha = 0.05, TVE=0.95,
                       distribution=c('norm'),
                       method = c('Welch','MC','Imhof'),
                       legend=TRUE,
                       highlight_changes=TRUE,
                       intervals = confidence_interval(x, CPs),
                       int.gradual=TRUE, ...
                       ){
  x <- funts(x)
  poss_types <- c('spaghetti','highdim', 'rainbow','banded','acf', 'pacf',
                  'summary', 'qq', 'distribution', 'change','interval',
                  'surface')
  type <- .verify_input(type, poss_types)

  if(length(CPs)==0) CPs <- NULL


  return_plot <- switch (type,
          spaghetti = {
            .plot_fd(X = x, CPs = CPs,
                                  plot_title = plot_title,
                                  val_axis_title = val_axis_title,
                                  res_axis_title = res_axis_title,
                                  FD_axis_title = FD_axis_title,
                                  eye = eye, aspectratio = aspectratio,
                                  showticklabels = showticklabels,
                                  interactive = TRUE)
          },
          highdim = {
            .plot_fd(X = x, CPs = CPs,
                                   plot_title = plot_title,
                                   val_axis_title = val_axis_title,
                                   res_axis_title = res_axis_title,
                                   FD_axis_title = FD_axis_title,
                                   eye = eye, aspectratio = aspectratio,
                                   showticklabels = showticklabels,
                                   interactive = FALSE)
          },
          rainbow = {
            .plot_rainbow(object = x, CPs=CPs)
          },
          banded = {
            .plot_banded(object=x, CPs=CPs, alpha=alpha)
          },
          acf = {
            # TODO:: Allow TVE
            acf.funts(x = x, lag.max = lag.max, alpha=alpha, method = method, ...)
          },
          pacf = {
            # TODO:: Allow TVE
            pacf.funts(x = x, lag.max = lag.max, alpha=alpha, n_pcs = NULL, method=method, ...)
          },
          summary = {
            summary.funts(object=x, CPs=CPs, lag.max=lag.max, d.max=d.max, ...)
          },
          qq = {
            qqplot(x=x, CPs=CPs, TVE=TVE, d.max=d.max,
                   alpha=alpha, legend=legend)
          },
          distribution = {
            .plot_distribution(X=x, CPs=CPs, TVE=TVE, d.max = d.max,
                     distribution = distribution,
                     alpha = alpha, legend = legend, ...)
          },
          change = {
            if(length(CPs)<1)
              stop('CPs must not be NULL to use change plot',call. = FALSE)
            .plot_change(X=x, CPs=CPs,
                        plot_title = plot_title,
                        val_axis_title = val_axis_title,
                        res_axis_title = res_axis_title,
                        FD_axis_title = FD_axis_title,
                        eye = eye,
                        aspectratio = aspectratio,
                        showticklabels = showticklabels,
                        warnings = FALSE)
          },
          interval = {
            .plot_interval(X=x, intervals=intervals,
                           plot_title = plot_title,
                           val_axis_title = val_axis_title,
                           res_axis_title = res_axis_title,
                           FD_axis_title = FD_axis_title,
                           eye = eye,
                           aspectratio = aspectratio,
                           showticklabels = showticklabels,
                           highlight_changes = highlight_changes,
                           int.gradual=int.gradual)
          },
          surface = {
            ## TODO:: Fix labels
            .plot_surface(x,...)
          },
          {
            stop('Review documentation for types of plots.',call. = FALSE)
          }
  )

  return_plot
}
