```{r plot_elec_paper,eval=FALSE}
plotData <- function(data){
  ### Define Variables
  curve_points = 1:nrow(data)
  plot_title = NULL
  val_axis_title = ""
  res_axis_title = ""
  FD_axis_title = ""
  FDReps = seq(as.Date("2014/01/01"), by = "day", length.out = 365)
  eye = list(x = -1.0, y = -3, z = 0.5)
  aspectratio=list(x=2,y=1.25,z=1.25)
  showticklabels=FALSE

  ## Plot
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

  tidyr::`%>%`(
    tidyr::`%>%`(
      tidyr::`%>%`(
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
              showticklabels = TRUE,
              dtick = "M2",
              tickformat="%b",
              tickfont =list(size=20)
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
}

plot_save <- plotData(electricity)

plotly::save_image(plot_save,
                   "C:/Users/jerem/Downloads/electricity_data.png",
                   width=1600, height=800)
```
