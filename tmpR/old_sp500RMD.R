

create_plot <- function(data){
  # https://plotly.com/r/time-series/
  ### Define Variables
  curve_points = 1:nrow(data)
  plot_title = NULL
  val_axis_title = ""
  res_axis_title = ""
  FD_axis_title = ""
  FDReps = colnames(data)
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
              dtick = "M6",
              tickformat="%b-%y",
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

create_plot1 <- function(data,CPs){
  curve_points = 1:nrow(data)
  plot_title = NULL
  val_axis_title = ""
  res_axis_title = ""
  FD_axis_title = ""
  FDReps = colnames(data)
  eye = list(x = -1.0, y = -3, z = 0.5)
  aspectratio=list(x=2,y=1.25,z=1.25)
  showticklabels=FALSE

  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- data.frame(
    "resolution" = NA,
    "FDRep" = NA,
    "Color" = NA,
    "Value" = NA
  )[-1, ]

  # Color Group to first CP
  for (j in 1:min(CPs)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[j],
        "Color" = 1,
        "Value" = data[, j]
      )
    )
  }
  # Color Group from last CP
  for (j in (max(CPs) + 1):number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = FDReps[j],
        "Color" = length(CPs) + 1,
        "Value" = data[, j]
      )
    )
  }

  # Color Additional Groups
  if (length(CPs) > 1) {
    for (i in 2:length(CPs)) {
      for (j in (CPs[i - 1] + 1):CPs[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = curve_points,
            "FDRep" = FDReps[j],
            "Color" = i,
            "Value" = data[, j]
          )
        )
      }
    }
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
  if (length(CPs) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
  }

  tidyr::`%>%`(
    tidyr::`%>%`(
      tidyr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                        type = "scatter3d", mode = "lines",
                        split = ~ as.factor(FDRep),
                        color = ~ as.factor(Color),
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
              dtick = "M6",
              tickformat="%b-%y",
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

SPYUS_ocidr <- compute_ocidr(SPYUS500)

set.seed(25641)
spy_min_cps <- binary_segmentation(SPYUS_ocidr,
                                   statistic='Tn',
                                   method='Sim',
                                   h_function = function(X){ncol(X)^(1/3)})
plot_ocidr <- create_plot(SPYUS_ocidr)
plotly::save_image(plot_save,
                   "C:/Users/jerem/Downloads/sp500_ocidr_min_1923.png",
                   width=1600, height=800)

tmp <- create_plot1(SPYUS_ocidr,spy_min_cps[[1]])
plotly::save_image(tmp,
                   "C:/Users/jerem/Downloads/sp500_ocidr_min_changes_1923.png",
                   width=1600, height=800)


<!-- ```{r} -->
  <!-- paperPlot <- function(data,CPs){ -->
      <!--   curve_points = 1:nrow(data) -->
        <!--   plot_title = NULL -->
          <!--   val_axis_title = "" -->
            <!--   res_axis_title = "" -->
              <!--   FD_axis_title = "" -->
                <!--   FDReps = colnames(data) -->
                  <!--   eye = list(x = -1.0, y = -3, z = 0.5) -->
                    <!--   aspectratio=list(x=2,y=1.25,z=1.25) -->
                      <!--   showticklabels=FALSE -->

                        <!--   number <- length(data[1, ]) -->
                          <!--   valRange <- c(min(data), max(data)) -->

                            <!--   plotData <- data.frame( -->
                                                             <!--     "resolution" = NA, -->
                                                             <!--     "FDRep" = NA, -->
                                                             <!--     "Color" = NA, -->
                                                             <!--     "Value" = NA -->
                                                             <!--   )[-1, ] -->

                              <!--   # Color Group to first CP -->
                              <!--   for (j in 1:min(CPs)) { -->
                                  <!--     plotData <- rbind( -->
                                                                <!--       plotData, -->
                                                                <!--       data.frame( -->
                                                                                         <!--         "resolution" = curve_points, -->
                                                                                         <!--         "FDRep" = FDReps[j], -->
                                                                                         <!--         "Color" = 1, -->
                                                                                         <!--         "Value" = data[, j] -->
                                                                                         <!--       ) -->
                                                                <!--     ) -->
                                    <!--   } -->
                              <!--   # Color Group from last CP -->
                              <!--   for (j in (max(CPs) + 1):number) { -->
                                  <!--     plotData <- rbind( -->
                                                                <!--       plotData, -->
                                                                <!--       data.frame( -->
                                                                                         <!--         "resolution" = curve_points, -->
                                                                                         <!--         "FDRep" = FDReps[j], -->
                                                                                         <!--         "Color" = length(CPs) + 1, -->
                                                                                         <!--         "Value" = data[, j] -->
                                                                                         <!--       ) -->
                                                                <!--     ) -->
                                    <!--   } -->

                              <!--   # Color Additional Groups -->
                              <!--   if (length(CPs) > 1) { -->
                                  <!--     for (i in 2:length(CPs)) { -->
                                      <!--       for (j in (CPs[i - 1] + 1):CPs[i]) { -->
                                          <!--         plotData <- rbind( -->
                                                                            <!--           plotData, -->
                                                                            <!--           data.frame( -->
                                                                                                         <!--             "resolution" = curve_points, -->
                                                                                                         <!--             "FDRep" = FDReps[j], -->
                                                                                                         <!--             "Color" = i, -->
                                                                                                         <!--             "Value" = data[, j] -->
                                                                                                         <!--           ) -->
                                                                            <!--         ) -->
                                            <!--       } -->
                                      <!--     } -->
                                  <!--   } -->

                              <!--   scene <- list( -->
                                                      <!--     camera = list(eye = eye), -->
                                                      <!--     aspectmode = "manual", -->
                                                      <!--     aspectratio = aspectratio -->
                                                      <!--   ) -->

                                <!--   # Get Colors -->
                                <!--   tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1") -->
                                  <!--   if (length(CPs) > 9) { -->
                                      <!--     tmpColors <- rep(tmpColors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)] -->
                                        <!--   } -->

                                  <!--   magrittr::`%>%`( -->
                                                            <!--     magrittr::`%>%`( -->
                                                                                        <!--       magrittr::`%>%`( -->
                                                                                                                      <!--         plotly::plot_ly(plotData, -->
                                                                                                                                                     <!--                         x = ~ as.factor(FDRep), y = ~resolution, z = ~Value, -->
                                                                                                                                                     <!--                         type = "scatter3d", mode = "lines", -->
                                                                                                                                                     <!--                         split = ~ as.factor(FDRep), -->
                                                                                                                                                     <!--                         color = ~ as.factor(Color), -->
                                                                                                                                                     <!--                         colors = tmpColors -->
                                                                                                                                                     <!--         ), -->
                                                                                                                      <!--         plotly::layout( -->
                                                                                                                                                     <!--           scene = list( -->
                                                                                                                                                                                    <!--             yaxis = list( -->
                                                                                                                                                                                                                     <!--               title = res_axis_title, -->
                                                                                                                                                                                                                     <!--               showticklabels = showticklabels -->
                                                                                                                                                                                                                     <!--             ), -->
                                                                                                                                                                                    <!--             xaxis = list( -->
                                                                                                                                                                                                                     <!--               title = FD_axis_title, -->
                                                                                                                                                                                                                     <!--               showticklabels = TRUE, -->
                                                                                                                                                                                                                     <!--               dtick = "M6",  -->
                                                                                                                                                                                                                     <!--               tickformat="%b-%y", -->
                                                                                                                                                                                                                     <!--               tickfont =list(size=20) -->
                                                                                                                                                                                                                     <!--             ), -->
                                                                                                                                                                                    <!--             zaxis = list( -->
                                                                                                                                                                                                                     <!--               title = val_axis_title, -->
                                                                                                                                                                                                                     <!--               showticklabels = showticklabels -->
                                                                                                                                                                                                                     <!--             ) -->
                                                                                                                                                                                    <!--           ) -->
                                                                                                                                                     <!--         ) -->
                                                                                                                      <!--       ), -->
                                                                                        <!--       plotly::layout(title = plot_title, scene = scene) -->
                                                                                        <!--     ), -->
                                                            <!--     plotly::layout(showlegend = FALSE) -->
                                                            <!--   ) -->
                                  <!-- } -->

  <!-- ``` -->
