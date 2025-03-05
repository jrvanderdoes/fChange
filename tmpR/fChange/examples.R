
X <- electricity

X$labels <- format(as.Date(X$labels),format='%b')

paper_plot <- function(X){
  val_axis_title = ''
  res_axis_title = ''
  FD_axis_title = ''
  plot_title = ''
  eye = list(x = -0.6, y = -2.5, z = 1.5)
  aspectratio = list(x = 1.75, y = 1, z = 1)

  plotData <- data.frame()
  for (i in 1:ncol(X)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$intratime,
        "FDRep" = i, #X$labels[i],
        "Value" = X$data[,i]
      )
    )
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
  plot_colors[6] <- "yellow"

  plt <- plotly::plot_ly(plotData,
                         x = ~FDRep, y = ~resolution, z = ~Value,
                         type = "scatter3d", mode = "lines",
                         color = ~ as.factor(FDRep),
                         colors = plot_colors
  ) %>%
    plotly::layout(
      scene = list(
        yaxis = list(
          title = res_axis_title,
          showticklabels = F,
          ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
          tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
        ),
        xaxis = list(
          title = FD_axis_title,
          showticklabels = T,
          ticktext=.select_n(vals=X$labels, n=5),
          tickvals=.select_n(vals=1:length(X$labels), n=5),
          tickfont = list(size = 15)
        ),
        zaxis = list(
          title = val_axis_title,
          showticklabels = F
        )
      )
    ) %>%
    plotly::layout(title = plot_title, scene = scene) %>%
    plotly::layout(showlegend = FALSE)
}


plotly::save_image(plt,'C:/Users/jerem/Downloads/data-can/electricity.png')
plotly::save_image(plt,'C:/Users/jerem/Downloads/data-can/cancer.png')

###########################################

X <- dfts(cancer)
val_axis_title = ''
res_axis_title = ''
FD_axis_title = ''
plot_title = ''
eye = list(x = -0.6, y = -2.5, z = 1.5)
aspectratio = list(x = 1.75, y = 1, z = 1)
# eye = list(x = -1, y = -1.5, z = 1.5)
# aspectratio = list(x = 1, y = 1, z = 1)

plotData <- data.frame()
for (i in 1:ncol(X)) {
  plotData <- rbind(
    plotData,
    data.frame(
      "resolution" = X$intratime,
      "FDRep" = i, #X$labels[i],
      "Value" = X$data[,i]
    )
  )
}

scene <- list(
  camera = list(eye = eye),
  aspectmode = "manual",
  aspectratio = aspectratio
)

plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
plot_colors[6] <- "yellow"

plt <- plotly::plot_ly(plotData,
                       x = ~FDRep, y = ~resolution, z = ~Value,
                       type = "scatter3d", mode = "lines",
                       color = ~ as.factor(FDRep),
                       colors = plot_colors
) %>%
  plotly::layout(
    scene = list(
      yaxis = list(
        title = res_axis_title,
        showticklabels = F,
        ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
        tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
      ),
      xaxis = list(
        title = FD_axis_title,
        showticklabels = T,
        ticktext=.select_n(vals=X$labels, n=5),
        tickvals=.select_n(vals=1:length(X$labels), n=5),
        tickfont = list(size = 15)
      ),
      zaxis = list(
        title = val_axis_title,
        showticklabels = F
      )
    )
  ) %>%
  plotly::layout(title = plot_title, scene = scene) %>%
  plotly::layout(showlegend = FALSE)

plotly::save_image(plt,'C:/Users/jerem/Downloads/data-can/cancer.png')

###########################################

X <- dfts(rates)
X$labels <- lubridate::year(as.Date(X$labels))

paper_plot1 <- function(X){
  changes <- NULL
  val_axis_title = ''
  res_axis_title = ''
  FD_axis_title = ''
  plot_title = ''
  aspectratio = c(2.5, .75, 1)


  # data <- X$data
  # curve_points <- X$intratime
  valRange <- c(
    floor(min(X$data,na.rm = T)),
    ceiling(max(X$data,na.rm = T))
  )

  name <- V1 <- value <- NULL
  data1 <- X$data
  colnames(data1) <- 1:ncol(X)
  plotData <- cbind(X$intratime,data1) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = 1+1:ncol(X)) %>%
    dplyr::mutate(name=as.numeric(name)) %>%
    dplyr::rename(resolution=V1,
                  Value=value,
                  FDRep=name)
  plotData <- plotData[order(plotData$FDRep),]

  ## Setup up Colors
  #   Rainbow for no changes, colored for changes
  plotData[["color"]] <- rep(1:ncol(X), each = nrow(X))
  if (!is.null(changes)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
    if (length(changes) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
    }

    changes <- unique(c(1, changes, ncol(X)))
    colors_plot <- rep(tmp_colors[1], ncol(X))
    for (i in 2:(length(changes) - 1)) {
      colors_plot[changes[i]:changes[i + 1]] <- tmp_colors[i]
    }
  } else {
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(X))
  }
  plotData[["color"]] <- rep(colors_plot, each = nrow(X) )

  ## Set up Tick Labels
  #   TODO:: Test more thoroughly
  r_z <- range(plotData$Value,na.rm = T)
  r_z[2] <- r_z[2]+4
  z_range <- round(r_z, -nchar(round(max(plotData$Value,na.rm = TRUE)))+1)

  scale_info <- list(
    col = "black", arrows = FALSE, cex = 1.5,
    x = list(
      at = seq(min(X$intratime),
               max(X$intratime),
               length.out = 5
      ),
      # labels=.specify_decimal(seq(max(X$intratime),
      #                                 min(X$intratime),
      #                                 length.out=5),2)
      labels = NULL
    ),
    y = list(
      at = -seq(1, ncol(X), length.out = 7),
      labels = X$labels[seq(1, ncol(X), length.out = 7)]
    ),
    z = list(
      at = seq(z_range[1], z_range[2],
               length.out = 5
      ),
      labels = NULL
    )
  )

  color <- NULL # Check note removal
  # lattice:::wireframe
  lattice::cloud(
    x = Value ~ (resolution) * (-FDRep),
    data = plotData,
    type = "l", groups = plotData$FDRep,
    par.box = c(col = "transparent"),
    par.settings =
      list(
        axis.line = list(col = "transparent"),
        superpose.line = list(col = (colors_plot))
      ),
    # screen=list(z = 90, x = -75,y=-45),
    # trellis.par.set(list(axis.text=list(cex=2)),
    #                "axis.line",list(col=NA)),
    xlim = rev(range(X$intratime)),
    zlim = valRange,
    aspect = aspectratio,
    drape = TRUE, colorkey = FALSE,
    scales = scale_info,
    xlab = list(res_axis_title, rot = 30,cex=1.75),
    ylab = list(paste0("\n", FD_axis_title), rot = -30,cex=1.75),
    zlab = list(val_axis_title, rot = 90, just = 0.75,cex=1.75),
    main = plot_title
  )
}


png('C:/Users/jerem/Downloads/data-can/rates.png')
plt
dev.off()

# ## Set up Tick Labels
# #   TODO:: Test more thoroughly
# r_z <- range(plotData$Value,na.rm = T)
# r_z[2] <- r_z[2]+4
# z_range <- round(r_z, -nchar(round(max(plotData$Value,na.rm = TRUE)))+1)
#
# scale_info <- list(
#   col = "black", arrows = FALSE, cex = 1.5,
#   x = list(
#     at = seq(min(X$intratime),
#              max(X$intratime),
#              length.out = 5
#     ),
#     # labels=.specify_decimal(seq(max(X$intratime),
#     #                                 min(X$intratime),
#     #                                 length.out=5),2)
#     labels = NULL
#   ),
#   y = list(
#     at = -seq(1, ncol(X), length.out = 7),
#     labels = X$labels[seq(1, ncol(X), length.out = 7)],
#     rot=45
#   ),
#   z = list(
#     at = seq(z_range[1], z_range[2],
#              length.out = 5
#     ),
#     labels = NULL
#   )
# )
#
# color <- NULL # Check note removal
# # lattice:::wireframe
# plt <- lattice::cloud(
#   x = Value ~ (resolution) * (-FDRep),
#   data = plotData,
#   type = "l", groups = plotData$FDRep,
#   par.box = c(col = "transparent"),
#   par.settings =
#     list(
#       axis.line = list(col = "transparent"),
#       superpose.line = list(col = (colors_plot))
#     ),
#   # screen=list(z = 90, x = -75,y=-45),
#   screen=list(z = 75, x = -60,y=-20),
#   # trellis.par.set(list(axis.text=list(cex=2)),
#   #                "axis.line",list(col=NA)),
#   xlim = rev(range(X$intratime)),
#   zlim = valRange,
#   aspect = aspectratio,
#   drape = TRUE, colorkey = FALSE,
#   scales = scale_info,
#   xlab = list(res_axis_title, rot = 60,cex=1.75),
#   ylab = list(paste0("\n", FD_axis_title), rot = -30,cex=1.75),
#   zlab = list(val_axis_title, rot = 90, just = 0.75,cex=1.75),
#   main = plot_title
# )
# plt
