#' Surface Plot
#'
#' Obtain a 3D plot of the data. This visualization is useful to detect any
#'  kind of dependency between the discretization points of the series.
#'
#' @param X Data
#' @param ... Arguments to pass to plotting. Includes: xlab, ylab, zlab, main,
#'  expand, ticktype, shade, col, theta and phi.
#'
#' @return Plot of the data
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_surface(autocovariance(electricity,1))
#'  }
.plot_surface <- function(X, ...){
  z <- t(X$data)

  # Color palette
  col.pal <- grDevices::colorRampPalette(c("blue", "red"))
  colors <- col.pal(100)

  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range <- cut(z.facet.center, 100)

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))     arguments$xlab <- ""
  if(!"ylab" %in% names(arguments))     arguments$ylab <- ""
  if(!"zlab" %in% names(arguments))     arguments$zlab <- ""
  if(!"main" %in% names(arguments))     arguments$main  <- X$name
  if(!"expand" %in% names(arguments))   arguments$expand  <- 0.6
  if(!"ticktype" %in% names(arguments)) arguments$ticktype  <-'detailed'
  if(!"shade" %in% names(arguments))    arguments$shade  <- NA
  if(!"col" %in% names(arguments))      arguments$col <- colors[z.facet.range]
  if(!"theta" %in% names(arguments))    arguments$theta  <- 25#315
  if(!"phi" %in% names(arguments))      arguments$phi  <- 30
  arguments$z <- z

  do.call(graphics::persp, arguments)
}
