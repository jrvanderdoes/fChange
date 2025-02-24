
number <- length(data[1, ])
valRange <- c(min(data), max(data))

plotData <- data.frame(
  "resolution" = NA,
  "FDRep" = NA,
  "Color" = NA,
  "Value" = NA
)[-1, ]

# Color Group to first CP
for (j in 1:min(changes)) {
  plotData <- rbind(
    plotData,
    data.frame(
      "resolution" = curve_points,
      "FDRep" = j,
      "Color" = 1,
      "Value" = data[, j]
    )
  )
}
# Color Group from last CP
for (j in (max(changes) + 1):number) {
  plotData <- rbind(
    plotData,
    data.frame(
      "resolution" = curve_points,
      "FDRep" = j,
      "Color" = length(changes) + 1,
      "Value" = data[, j]
    )
  )
}

# Color Additional Groups
if (length(changes) > 1) {
  for (i in 2:length(changes)) {
    for (j in (changes[i - 1] + 1):changes[i]) {
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = curve_points,
          "FDRep" = j,
          "Color" = i,
          "Value" = data[, j]
        )
      )
    }
  }
}

# Setup Changes / Confidence Intervals
ci_data_new <- ci_data
ci_data_new$lower <- floor(ci_data_new$lower)
ci_data_new$upper <- ceiling(ci_data_new$upper)

res <- length(unique(plotData$resolution))
# Color Additional Groups
for (i in 1:nrow(ci_data_new)) {
  change <- ci_data_new$change[i]
  prev_color <- unique(plotData[plotData$FDRep==(change-1),"Color"])
  next_color <- unique(plotData[plotData$FDRep==(change+1),"Color"])
  plotData[plotData$FDRep==change,"Color"] <- '#000000'
  # Lower Ramp
  for (j in (change-1):ci_data_new[i,"lower"]) {
    # curr_color <- plotData[plotData$FDRep==j & plotData$resolution>res/2,"Color"]
    plotData[plotData$FDRep==j,"Color"] <- paste0(prev_color,'.25')
  }
  # Upper Ramp
  for (j in (change+1):ci_data_new[i,"upper"]) {
    # curr_color <- plotData[plotData$FDRep==j,"Color"]
    plotData[plotData$FDRep==j,"Color"] <- paste0(next_color,'.75')
  }
}

scene <- list(
  camera = list(eye = eye),
  aspectmode = "manual",
  aspectratio = aspectratio
)

# Get Colors
tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
if (length(changes) > 9) {
  tmpColors <- rep(tmpColors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
}

finalcolors <- c(tmpColors[1], colorspace::lighten(tmpColors[1],amount=0.5))
for(i in 2:length(tmpColors)) {
  finalcolors <- c(finalcolors, colorspace::lighten(tmpColors[i],amount=0.5),
                   tmpColors[i], colorspace::lighten(tmpColors[i],amount=0.5))
}
tmpColors <- c('black',tmpColors)

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
