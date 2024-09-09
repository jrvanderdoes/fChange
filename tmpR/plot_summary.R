#' Title
#'
#' @param X
#'
#' @return
#' @export
#'
#' @examples
#' .summary_nochange(
#'   generate_brownian_motion(N = 500,
#'      v=seq(from = 0, to = 1, length.out = 20)))
#' .summary_nochange(funts(electricity))
#' .summary_nochange(funts(generateFAR1(20,300)))
.summary_nochange <- function(X){

  #####
  ## Plot white noise lags
  wn_pvalues <- rep(NA, max.lag)
  for(h in 1:max.lag){
    wn_pvalues[h] <- single_lag_test(X,lag = h,method = 'bootstrap')$p_value
  }

  plot_whitenoise <-
    ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=1:max.lag, y=wn_pvalues)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),
                          linetype='dotted', col='red') +
      ggplot2::theme_bw() +
      ggplot2::ylim(c(0,1)) +
      ggplot2::xlab('') +
      ggplot2::ylab('')

  #####
  ## Plot rainbow TODO:: QQ-plot??
  plot_rainbow <- .plot_stack(X)

  #####
  ## PlotA ACF
  data_acf <- acf(X,figure=FALSE)

  plot_acf <-
    ggplot2::ggplot(
      mapping=ggplot2::aes( x=1:length(data_acf$WWN_bound) )
    ) +
    ggplot2::geom_segment( ggplot2::aes(y=0, yend=data_acf$rho)) +
    ggplot2::geom_line(ggplot2::aes( y=data_acf$WWN_bound),
                       col='red', linetype='dashed'  ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=data_acf$SWN_bound),
                        col='blue', linetype='longdash' ) +
    ggplot2::theme_bw() +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot Lines
  data_lines <- cbind(data.frame('time'=X$intraobs),
               X$data) %>%
    pivot_longer(cols = 1+1:ncol(X$data))

  colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
  colors_plot[6] <- "yellow"
  colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(X$data))
  data_lines$color <- rep(colors_plot, nrow(X$data) )

  # colors <- as.character(1:ncol(X$data))
  # for(i in 1:length(colors)){
  #   if(nchar(colors[i]) < max(nchar(colors)) ){
  #     colors[i] <- paste0(
  #       paste0(rep('0', max(nchar(colors))-nchar(colors[i])),collapse=""),
  #       colors[i])
  #   }
  # }
  # data_lines$color <- rep(colors,times=nrow(X$data))

  # colors <- as.character(1:ncol(X$data))
  # for(i in 1:length(colors)){
  #   if(nchar(colors[i]) < max(nchar(colors)) ){
  #     colors[i] <- paste0(
  #       paste0(rep('0', max(nchar(colors))-nchar(colors[i])),collapse=""),
  #       colors[i])
  #   }
  # }
  # data_lines$color <- rep(colors,times=nrow(X$data))
  data_lines$name <- as.numeric(data_lines$name)

  plot_lines <- ggplot2::ggplot(data_lines,
                          aes(y=time,
                              x=name,#as.Date(name)),
                              z=value,
                              color=color)) +#as.character(color))) +
    theme_void() +
    stat_3D(theta=0, phi=15, geom='path') +
    ggplot2::scale_color_manual(
      breaks = data_lines$color,
      values = data_lines$color
    ) +
    # ggplot2::scale_fill_gradientn(
    #   limits  = range(0,max(data_lines$value)),
    #   colours = c('red','yellow','blue')) +
    ggplot2::guides(color='none')

  #####
  ## Plot Summary
  ggpubr::ggarrange(plot_lines,
                    ggpubr::ggarrange(plot_acf,
                                      plot_rainbow,
                                      ncol = 2),
                    plot_whitenoise,
                    nrow = 3 )
}


#' Title
#'
#' @param X
#' @param CPs
#'
#' @return
#' @export
#'
#' @examples
.plot_stack <- function(X, CPs=NULL){
  data <- X$data
  plot_data <- tidyr::pivot_longer(data = cbind(data.frame('Time'=X$intraobs),data),
                      cols = 1+1:ncol(data))
  # Colors
  ## Setup up Colors
  #   Rainbow for no CPs, colored for CPs
  # col_data <- unique(plot_data[,"name"])
  # col_data$color <- 1:nrow(col_data)
  # plot_data[["color"]] <- rep(1:ncol(data), times = nrow(data))
  if (!is.null(CPs)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
    if (length(CPs) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
    }

    CPs <- unique(c(1, CPs, ncol(data)))
    colors_plot <- rep(tmp_colors[1], ncol(data))
    for (i in 2:(length(CPs) - 1)) {
      colors_plot[CPs[i]:CPs[i + 1]] <- tmp_colors[i]
    }
  } else {
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(data))
    # plot_data[["color"]] <- rep(colors_plot,
    #                            times = nrow(data)
    #)
  }
  # col_data$color <- colors_plot
  plot_data$color <- rep(colors_plot,times=nrow(data))

  # Order
  max_val <- max(nchar(plot_data$name))
  for(i in 1:nrow(plot_data)){
    if(nchar(plot_data$name[i]) < max_val ){
      plot_data$name[i] <- paste0(
        paste0(rep('0', max_val-nchar(plot_data$name[i])),collapse=""),
        plot_data$name[i])
    }
  }

  # Plot
  ggplot2::ggplot(plot_data) +
    ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
                                    group=name,
                                    color=color),
                                    # alpha=0.5 + 0.1 *
                                    #   (1 - as.numeric(as.factor(color)) /
                                    #      max(as.numeric(as.factor(color))))),
                                    # alpha=0.75),#alpha=0.5,
                                    # alpha=0.95 + 0.05 *
                                    #   (1-as.numeric(name)/max(as.numeric(name)))),
                       alpha=0.5, linetype='dashed') +
    ggplot2::theme_bw() +
    ggplot2::guides(color="none",alpha='none') +
    ggplot2::scale_color_manual(values = plot_data$color,
                                breaks = plot_data$color) +
    ggplot2::xlab('') +
    ggplot2::ylab('')
    # ggplot2::scale_colour_gradientn(name = 'category', colours = c('blue', 'red'))
    #ggplot2::scale_colour_manual(values=cols)
}


# plot_single_lags <- function(X, max.lag=20){
#   # Demean based on CP
#
#   p_values <- rep(NA, max.lag)
#   for(h in 1:max.lag){
#     p_values[h] <- single_lag_test(X,lag = h,method = 'bootstrap')$p_value
#   }
#
#   ggplot2::ggplot() +
#     ggplot2::geom_point(ggplot2::aes(x=1:max.lag, y=p_values)) +
#     ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),linetype='dotted', col='red') +
#     ggplot2::theme_bw() +
#     ggplot2::ylim(c(0,1)) +
#     ggplot2::xlab('') +
#     ggplot2::ylab('')
# }


#' Draw 3D Geoms
#'
#' This function adds 3D geoms such as points and paths to a ggplot2 plot.
#'
#' @param theta The azimuthal direction in degrees.
#' @param phi The colatitude in degrees.
#' @param geom The geom type to use *ie. "point", "path", "line"*
#' @param ... Arguements passed on to layer.
#' These are often aesthetics, used to set an
#' aesthetic to a fixed value, like color = "red" or size = 3.
#'
#' @references Acker D (2024). gg3D: 3D perspective plots for ggplot2. R
#'                package version 0.0.0.9000.
stat_3D <- function(mapping = NULL, data = NULL, geom = "point",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
  # 3d Plot prototype
  Stat3D <- ggplot2::ggproto("Stat3D",  ggplot2::Stat,
                             setup_params = function(data, params) {

                               params$xrange <- range(data$x)
                               params$yrange <- range(data$y)
                               params$zrange <- range(data$z)

                               params
                             },
                             compute_group = function(data, scales,
                                                      theta=135, phi=60,
                                                      xrange=c(0,1),
                                                      yrange=c(0,1),
                                                      zrange=c(0,1)) {
                               data <- data %>%
                                 dplyr::mutate(
                                   x = scales::rescale(x, from=xrange, to=c(0,1)),
                                   y = scales::rescale(y, from=yrange, to=c(0,1)),
                                   z = scales::rescale(z, from=zrange, to=c(0,1)))

                               pmat <- plot3D::perspbox(z=diag(2), plot=FALSE,
                                                       theta=theta, phi=phi)

                               XY <- plot3D::trans3D(
                                 x = data$x,
                                 y = data$y,
                                 z = data$z,
                                 pmat = pmat) %>%
                                 data.frame()

                               data$x <- XY$x
                               data$y <- XY$y

                               data
                             },
                             required_aes = c("x", "y", "z")
  )

  layer(
    stat = Stat3D, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
