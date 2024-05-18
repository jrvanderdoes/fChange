
library(plotly)
library(RColorBrewer)
library(magrittr)

#' Plot Functional Time Series
#'
#' This is a base function to plot functional time series. There are a lot of
#'  possible settings to change in the function. Happy to help if anything is
#'  amiss as I had to modify a similar function to format it properly for our
#'  paper.
#'
#' @param data Data.frame of data with each column being an observation
#' @param curve_points The points each observation is seen.
#'    Typically this is evenly spaced and thus the default is sufficient.
#' @param plot_title String for the title of the plot
#' @param val_axis_title String for the axis title of values
#' @param res_axis_title String for the axis title of resolution
#' @param FD_axis_title String for the axis title of observations
#' @param FDReps Vector of observation labels.
#'    You can convert this to dates as well, but you may need to update the function
#' @param eye List giving the angle the plot is viewed
#' @param aspectratio List giving the ratio for each axis
#' @param showticklabels Boolean to hide some tick labels
#' @param save_path String that gives the path to export the figure.
#'    This is much higher resolution than exporting from Rstudio.
#'    Python and several packages are required (like plotly and perhaps others)
#'
#'
#' @return
#' @export
#'
#' @examples
plot_fd <- function(data,
                    curve_points = 1:nrow(data),
                    plot_title = NULL,
                    val_axis_title = "Value",
                    res_axis_title = "Resolution",
                    FD_axis_title = "Observation",
                    FDReps = 1:ncol(data),
                    eye = list(x = -1.5, y = -1.5, z = 1.5),
                    aspectratio = list(x = 1, y = 1, z = 1),
                    showticklabels = TRUE,
                    save_path = NULL) {
  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- data.frame(
    "resolution" = curve_points,
    "FDRep" = FDReps[1],
    "Value" = data[, 1]
  )

  for (i in 2:number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[i],
        "Value" = as.factor(data[, i])
      )
    )
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  tmpColors <- RColorBrewer::brewer.pal(11, "Spectral")
  tmpColors[6] <- "yellow"

  export_plot <- magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~FDRep, y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        color = ~ as.factor(FDRep),
                        colors = tmpColors
        ),
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
        )
      ),
      plotly::layout(title = plot_title, scene = scene)
    ),
    plotly::layout(showlegend = FALSE)
  )

  if(!is.null(save_path)){
    plotly::save_image(export_plot,
                       save_path,
                       width=800, height=800)
  }

  export_plot
  # plotly::plot_ly(plotData,
  #                x = ~as.factor(FDRep), y = ~resolution, z = ~Value,
  #                type = 'scatter3d', mode = 'lines',
  #                color = ~as.factor(FDRep),
  #                colors = tmpColors) %>%
  #  #colors = c("red", "yellow", "blue")) %>%
  #  #colors='Spectral') %>%
  #  plotly::layout(
  #    scene = list(
  #      yaxis = list(title = "resolution"),
  #      xaxis = list(title = "FD Sims"),
  #      zaxis = list(title = "Value")
  #    )) %>%
  #  plotly::layout(title = plot_title, scene = scene) %>%
  #  plotly::layout(showlegend = FALSE)
  ## line = list(width = 4, color = ~as.factor(FDRep),
  ##            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
}
