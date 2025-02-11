
#' Plot Data with Forecast Bands
#'
#' @inheritParams plot.funts
#' @param lower Matrix of resolution-by-estimates indicating the lower bounds
#' @param upper Matrix of resolution-by-estimates indicating the upper bounds
#'
#' @returns Plotly Plot
#'
#' @noRd
#' @keywords internal
.plot_forecast <- function(X, lower, upper, CPs=NULL,
                           plot_title = X$name,
                           val_axis_title = "Value",
                           res_axis_title = "Resolution",
                           FD_axis_title = "Observations",
                           eye = list(x = -1.5, y = -1.5, z = 1.5),
                           aspectratio = NULL,
                           showticklabels = TRUE){
  # Get Sizes
  pred_n <- max(ncol(lower),ncol(upper))
  CPs <- unique(c(CPs, ncol(X)-pred_n))

  # Setup Data
  plotData <- plotData_forecast <- plotData_lower <- plotData_upper <- data.frame()

  # Main Data
  if (length(CPs) > 1) {
    for (j in 1:min(CPs)) {
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$intraobs,
          "FDRep" = j,
          "Color" = 1,
          "Value" = X$data[, j]
        )
      )
    }

    for (i in 2:length(CPs)) {
      for (j in (CPs[i - 1] + 1):CPs[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = X$intraobs,
            "FDRep" = j,
            "Color" = i,
            "Value" = X$data[, j]
          )
        )
      }
    }
  } else{
    for (j in 1:min(CPs)) {
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$intraobs,
          "FDRep" = j,
          "Color" = j,
          "Value" = X$data[, j]
        )
      )
    }

  }

  # Forecast Data
  for (j in (max(CPs) + 1):ncol(X)) {
    plotData_forecast <- rbind(
      plotData_forecast,
      data.frame(
        "resolution" = X$intraobs,
        "FDRep" = j,
        "Color" = max(CPs) + 1,
        "Value" = X$data[, j]
      )
    )

    plotData_lower <- rbind(
      plotData_lower,
      data.frame(
        "resolution" = X$intraobs,
        "FDRep" = j,
        "Color" = max(CPs) + 1,
        "Value" = lower[,j-max(CPs)]
      )
    )

    plotData_upper <- rbind(
      plotData_upper,
      data.frame(
        "resolution" = X$intraobs,
        "FDRep" = j,
        "Color" = max(CPs) + 1,
        "Value" = upper[,j-max(CPs)]
      )
    )
  }


  # Get Colors
  if (length(CPs) > 1) {
    plot_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs))), "Set1")[1:min(9,length(CPs))]
    if (length(CPs)-1 > 9) {
      plot_colors <- rep(plot_colors, ceiling(c(length(CPs) + 1) / 9))[1:length(CPs)]
    }

    plot_colors <- c(plot_colors,'black')
    names(plot_colors) <- c(1:(length(CPs)),CPs[length(CPs)]+1)
  } else{
    plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
    plot_colors[6] <- "yellow"
    plot_colors <- c(grDevices::colorRampPalette(plot_colors)(CPs),'black')
    names(plot_colors) <- c(1:(CPs+1))
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
