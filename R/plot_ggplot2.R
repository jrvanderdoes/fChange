#' Draw 3D Geoms for ggplot2
#'
#' This function adds 3D geoms such as points and paths to a ggplot2 plot.
#'
#' @inheritParams ggplot2::ggplot
#' @inheritParams ggplot2::layer
#' @param position Set position information. Find more details in \code{ggplot2}
#'  function.
#' @param na.rm Boolean if na data should be removed.
#' @param ... Arguments passed on to layer. Often the aesthetics
#'  like \code{color = "red"} or \code{size = 3}. Two important ones are
#'  \code{theta} (azimuthal rotation) and \code{phi} (colatitude rotation)
#'  to define angle in degrees of viewing data.
#'
#' @references Acker D (2024). gg3D: 3D perspective plots for ggplot2. R
#'                package version 0.0.0.9000.
#'
#' @export
#' @returns No direct return value, called to be used with [ggplot2::ggplot()]
#'  in designing the plot.
#'
#' @examples
#'dat <- electricity
#'
#'data_lines <- cbind(data.frame('Time'=dat$fparam), dat$data) %>%
#' tidyr::pivot_longer(cols = 1+1:ncol(dat$data))
#'
#'colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
#'colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(dat$data))
#'data_lines$color <- rep(colors_plot, nrow(dat$data) )
#'data_lines$name <- as.numeric(data_lines$name)
#'
#'result <- ggplot2::ggplot(data_lines,
#'   ggplot2::aes(y=Time, x=name, z=value, color=color)) +
#'   ggplot2::theme_void() +
#'   stat_3D(theta=0, phi=15, geom='path') +
#'   ggplot2::scale_color_manual(
#'      breaks = data_lines$color,
#'      values = data_lines$color
#'   ) +
#'   ggplot2::guides(color='none')
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
                               # Make NULL to remove warnings
                               x <- y <- z <- NULL
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

  ggplot2::layer(
    stat = Stat3D, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
