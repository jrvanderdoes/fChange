
#' Plot Data with Forecast Bands
#'
#' @inheritParams plot.dfts
#' @param lower Matrix of resolution-by-estimates indicating the lower bounds
#' @param upper Matrix of resolution-by-estimates indicating the upper bounds
#'
#' @returns Plotly Plot
#'
#' @noRd
#' @keywords internal
.plot_forecast <- function(X, lower, upper, changes=NULL,
                           plot_title = X$name,
                           val_axis_title = "Value",
                           res_axis_title = "Resolution",
                           FD_axis_title = "Observations",
                           eye = list(x = -1.5, y = -1.5, z = 1.5),
                           aspectratio = NULL,
                           showticklabels = TRUE){
  # Get Sizes
  pred_n <- max(ncol(lower),ncol(upper))
  changes <- unique(c(changes, ncol(X)-pred_n))

  # Setup Data
  plotData <- plotData_forecast <- plotData_lower <- plotData_upper <- data.frame()

  # Main Data
  if (length(changes) > 1) {
    for (j in 1:min(changes)) {
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$fparam,
          "FDRep" = j,
          "Color" = 1,
          "Value" = X$data[, j]
        )
      )
    }

    for (i in 2:length(changes)) {
      for (j in (changes[i - 1] + 1):changes[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = X$fparam,
            "FDRep" = j,
            "Color" = i,
            "Value" = X$data[, j]
          )
        )
      }
    }
  } else{
    for (j in 1:min(changes)) {
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$fparam,
          "FDRep" = j,
          "Color" = j,
          "Value" = X$data[, j]
        )
      )
    }

  }

  # Forecast Data
  for (j in (max(changes) + 1):ncol(X)) {
    plotData_forecast <- rbind(
      plotData_forecast,
      data.frame(
        "resolution" = X$fparam,
        "FDRep" = j,
        "Color" = max(changes) + 1,
        "Value" = X$data[, j]
      )
    )

    plotData_lower <- rbind(
      plotData_lower,
      data.frame(
        "resolution" = X$fparam,
        "FDRep" = j,
        "Color" = max(changes) + 1,
        "Value" = lower[,j-max(changes)]
      )
    )

    plotData_upper <- rbind(
      plotData_upper,
      data.frame(
        "resolution" = X$fparam,
        "FDRep" = j,
        "Color" = max(changes) + 1,
        "Value" = upper[,j-max(changes)]
      )
    )
  }


  # Get Colors
  if (length(changes) > 1) {
    plot_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes))), "Set1")[1:min(9,length(changes))]
    if (length(changes)-1 > 9) {
      plot_colors <- rep(plot_colors, ceiling(c(length(changes) + 1) / 9))[1:length(changes)]
    }

    plot_colors <- c(plot_colors,'black')
    names(plot_colors) <- c(1:(length(changes)),changes[length(changes)]+1)
  } else{
    plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
    plot_colors[6] <- "yellow"
    plot_colors <- c(grDevices::colorRampPalette(plot_colors)(changes),'black')
    names(plot_colors) <- c(1:(changes+1))
  }

  # Setup View
  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  # Plot
  plotly::plot_ly(rbind(plotData,plotData_forecast),
                  x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  split = ~ as.factor(FDRep),
                  color = ~ as.factor(Color),
                  colors = plot_colors
  ) %>%
    plotly::add_lines(data=plotData_lower,
                      x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                      type = "scatter3d", mode = "lines",
                      split = ~ as.factor(FDRep),
                      colors = 'black',
                      opacity=0.15
    ) %>%
    plotly::add_lines(data=plotData_upper,
                      x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                      type = "scatter3d", mode = "lines",
                      split = ~ as.factor(FDRep),
                      colors = 'black',
                      opacity=0.15
    ) %>%
    plotly::layout(
      scene = list(
        yaxis = list(
          title = res_axis_title,
          showticklabels = showticklabels
        ),
        xaxis = list(
          title = FD_axis_title,
          showticklabels = showticklabels
        ),
        zaxis = list(
          title = val_axis_title,
          showticklabels = showticklabels
        )
      )
    ) %>%
    plotly::layout(title = plot_title, scene = scene) %>%
    plotly::layout(showlegend = FALSE)
}
