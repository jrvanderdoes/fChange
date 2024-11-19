curve_points <- 1:nrow(data)

number <- ncol(data)
valRange <- c(min(data), max(data))

plotData <- data.frame(
  "resolution" = NA,
  "FDRep" = NA,
  "Color" = NA,
  "Value" = NA
)[-1, ]

# Color Group to first CP
for (j in 1:number) {
  plotData <- rbind(
    plotData,
    data.frame(
      "resolution" = curve_points,
      "FDRep" = j,
      "Color" = '#d3d3d3',
      "Value" = data[, j]
    )
  )
}

# Setup Changes / Confidence Intervals
ci_data_new <- ci_data
ci_data_new$lower <- floor(ci_data_new$lower)
ci_data_new$upper <- ceiling(ci_data_new$upper)

# Color Additional Groups
for (i in 1:nrow(ci_data_new)) {
  change <- ci_data_new$change[i]
  plotData[plotData$FDRep==change,"Color"] <- '#ff0000'
  # Lower Ramp
  color_ramp <-
    colorRampPalette(c('#de6fa1','#ffe4e1'))(#c('#208fbc','#d8eff8'))(
      change-ci_data_new[i,"lower"])
  for (j in (change-1):ci_data_new[i,"lower"]) {
    curr_color <- plotData[plotData$FDRep==j,"Color"]
    plotData[plotData$FDRep==j,"Color"] <-
      ifelse(curr_color!='#d3d3d3',
         min(curr_color, color_ramp[change-j]),
         color_ramp[change-j])
  }
  # Upper Ramp
  color_ramp <-
    colorRampPalette(c('#de6fa1','#ffe4e1'))(
      ci_data_new[i,"upper"] - change)
  for (j in (change+1):ci_data_new[i,"upper"]) {
    curr_color <- plotData[plotData$FDRep==j,"Color"]
    plotData[plotData$FDRep==j,"Color"] <-
      ifelse(curr_color!='#d3d3d3',
         min(curr_color, color_ramp[j-change]),
         color_ramp[j-change])
  }
}
scene <- list(
  camera = list(eye = eye),
  aspectmode = "manual",
  aspectratio = aspectratio
)

# # Get Colors
tmpColors <- levels(as.factor(plotData$Color))
# tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
# if (length(CPs) > 9) {
#   tmpColors <- rep(tmpColors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
# }

magrittr::`%>%`(
  magrittr::`%>%`(
    magrittr::`%>%`(
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
