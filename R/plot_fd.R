#' Plot functional data
#'
#' \code{.plot_fd} plots functional data either in an fd object or evaluated at certain points
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param changes (Optional) Vectors of numeric values indicating the location of
#'     change points. This will color each section differently. Default vallue
#'     is NULL.
#' @param plot_title (Optional) String to title the plot. Default value is NULL.
#' @param val_axis_title (Optional) String to title the axis for values of
#'     observation. Default value is 'Value'.
#' @param res_axis_title (Optional) String to title the axis for the resolution
#'     range axis (observations in an FD object). Default value is 'resolution'.
#' @param FD_axis_title (Optional) String to title the axis for FD observations.
#'     Default value is 'Observation'.
#' @param eye (Optional) List with certain parameters to determine the view of
#'     the resulting image. The default value is list(x = -1.5, y = -1.5, z = 1.5).
#' @param aspectratio (Optional) List with certain parameters to determine the
#'     image size / viewing parameters. The default value is
#'     list(x=1, y=1, z=1), which is used in the non-interactive case. If
#'     interactive is TRUE, then the default is c(2.5,.75,1) and all lists are
#'     converted to this value so enter a vector if you wish to modify it for
#'     the non-interactive case.
#' @param showticklabels (Optional) Boolean to indicate if the tick labels should
#'     be given in the image. Default value is TRUE.
#' @param interactive (Optional) Boolean value to indicate if the plot should
#'  be interactive or not. Recommended to be TRUE when the number of observations
#'  is high. Default is FALSE
#'
#' @return A plot for the data. It is an interactive plotly plot if
#'  interactive is FALSE (default) and a lattice plot if TRUE.
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_fd(X = electricity$data[, 1:10],interactive = FALSE)
#'    \item .plot_fd(X = electricity$data[, 1:50], changes = c(25))
#'  }
.plot_fd <- function(X, changes = NULL, plot_title = X$name,
                    val_axis_title = "Value",
                    res_axis_title = "fparam",
                    FD_axis_title = "Observations",
                    eye = list(x = -1.5, y = -1.5, z = 1.5),
                    aspectratio = NULL,
                    showticklabels = TRUE, interactive = TRUE) {
  # Setup
  if(is.null(eye)) eye <- list(x = -1.5, y = -1.5, z = 1.5)
  if(is.null(val_axis_title)) val_axis_title = ''
  if(is.null(res_axis_title)) res_axis_title = ''
  if(is.null(FD_axis_title)) FD_axis_title = ''

  X <- dfts(X,inc.warnings = FALSE)
  if(!is.null(changes)) changes <- changes[order(changes)]

  if (!interactive) {
    if(is.null(aspectratio)) aspectratio <- c(2.5, 0.75, 1)

    fdPlot <- .plot_evalfd_highdim(
      X = X,
      changes = changes, plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title,
      aspectratio = aspectratio,
      showticklabels = showticklabels
    )
  } else if (!is.null(changes) && length(stats::na.omit(changes)) > 0) {
    if(is.null(aspectratio)) aspectratio <- list(x = 1, y = 1, z = 1)

    fdPlot <- .plot_evalfd_3dlines_changes(
      X = X,
      changes = changes, plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title,
      eye = eye, aspectratio = aspectratio,
      showticklabels = showticklabels
    ) # ,FDReps)
  } else {
    if(is.null(aspectratio)) aspectratio <- list(x = 1, y = 1, z = 1)

    fdPlot <- .plot_evalfd_3dlines(
      X=X,
      plot_title = plot_title,
      val_axis_title = val_axis_title,
      res_axis_title = res_axis_title,
      FD_axis_title = FD_axis_title,
      eye = eye, aspectratio = aspectratio,
      showticklabels = showticklabels
    )
  }

  fdPlot
}


#' Plot Functional Data Without Change Points
#'
#' This (internal) function to plot the function data, with no coloring based on
#'     change points.
#'
#' @inheritParams .plot_fd
#'
#' @return A plotly plot
#'
#' @noRd
#' @keywords internal
.plot_evalfd_3dlines <- function(X, plot_title = NULL,
                                 val_axis_title = "Value",
                                 res_axis_title = "fparam",
                                 FD_axis_title = "Observations",
                                 eye = list(x = -1.5, y = -1.5, z = 1.5),
                                 aspectratio = list(x = 1, y = 1, z = 1),
                                 showticklabels = TRUE) {

  plotData <- data.frame()
  for (i in 1:ncol(X)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$fparam,
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

  plotly::plot_ly(plotData,
                  x = ~FDRep, y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  color = ~ as.factor(FDRep),
                  colors = plot_colors
  ) %>%
    plotly::layout(
      scene = list(
        yaxis = list(
          title = res_axis_title,
          showticklabels = showticklabels,
          ticktext=round(seq(min(X$fparam),max(X$fparam),length.out=5),3),
          tickvals=.select_n(vals=seq(min(X$fparam),max(X$fparam),length.out=5), n=5)
        ),
        xaxis = list(
          title = FD_axis_title,
          showticklabels = showticklabels,
          ticktext=.select_n(vals=X$labels, n=6),
          tickvals=.select_n(vals=1:length(X$labels), n=6)
        ),
        zaxis = list(
          title = val_axis_title,
          showticklabels = showticklabels
        )
      )
    ) %>%
    plotly::layout(title = plot_title, scene = scene) %>%
    plotly::layout(showlegend = FALSE)

  # magrittr::`%>%`(
  #   magrittr::`%>%`(
  #     magrittr::`%>%`(
  #       plotly::plot_ly(plotData,
  #         x = ~FDRep, y = ~resolution, z = ~Value,
  #         type = "scatter3d", mode = "lines",
  #         color = ~ as.factor(FDRep),
  #         colors = plot_colors
  #       ),
  #       plotly::layout(
  #         scene = list(
  #           yaxis = list(
  #             title = res_axis_title,
  #             showticklabels = showticklabels
  #           ),
  #           xaxis = list(
  #             title = FD_axis_title,
  #             showticklabels = showticklabels
  #           ),
  #           zaxis = list(
  #             title = val_axis_title,
  #             showticklabels = showticklabels
  #           )
  #         )
  #       )
  #     ),
  #     plotly::layout(title = plot_title, scene = scene)
  #   ),
  #   plotly::layout(showlegend = FALSE)
  # )

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


#' Plot Functional Data With Change Points
#'
#' This (internal) function to plot the function data, with coloring based on
#'     change points.
#'
#' @inheritParams .plot_fd
#' @param changes Vectors of numeric values indicating the location of change points.
#'     This will color each section differently.
#'
#' @return A plotly plot
#'
#' @noRd
#' @keywords internal
.plot_evalfd_3dlines_changes <- function(X, changes,
                                     plot_title = NULL,
                                     val_axis_title = "Value",
                                     res_axis_title = "fparam",
                                     FD_axis_title = "Observations",
                                     eye = list(x = -1.5, y = -1.5, z = 1.5),
                                     aspectratio = list(x = 1, y = 1, z = 1),
                                     showticklabels = TRUE) {
  plotData <- data.frame()

  # Color Group to first CP
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
  # Color Group from last CP
  for (j in (max(changes) + 1):ncol(X)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$fparam,
        "FDRep" = j,
        "Color" = length(changes) + 1,
        "Value" = X$data[, j]
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
            "resolution" = X$fparam,
            "FDRep" = j,
            "Color" = i,
            "Value" = X$data[, j]
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
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
  if (length(changes) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
  }

  # magrittr::`%>%`(
  #   magrittr::`%>%`(
  #     magrittr::`%>%`(
  #       plotly::plot_ly(plotData,
  #         x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
  #         type = "scatter3d", mode = "lines",
  #         split = ~ as.factor(FDRep),
  #         color = ~ as.factor(Color),
  #         colors = tmpColors
  #       ),
  #       plotly::layout(
  #         scene = list(
  #           yaxis = list(
  #             title = res_axis_title,
  #             showticklabels = showticklabels
  #           ),
  #           xaxis = list(
  #             title = FD_axis_title,
  #             showticklabels = showticklabels
  #           ),
  #           zaxis = list(
  #             title = val_axis_title,
  #             showticklabels = showticklabels
  #           )
  #         )
  #       )
  #     ),
  #     plotly::layout(title = plot_title, scene = scene)
  #   ),
  #   plotly::layout(showlegend = FALSE)
  # )
  plotly::plot_ly(plotData,
                  x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  split = ~ as.factor(FDRep),
                  color = ~ as.factor(Color),
                  colors = tmpColors
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


#' Plot With Surface for Speed
#'
#' @inheritParams .plot_fd
#' @param aspectratio (Optional) List with certain parameters to determine the
#'     image size parameters. The default value is c(2.5,.75,1). Also any lists
#'     are converted to this (as it is likely from default call in
#'     \code{.plot_fd()}).
#'
#' @return A lattice plot
#'
#' @noRd
#' @keywords internal
.plot_evalfd_highdim <- function(X, changes = NULL,
                                 plot_title = NULL,
                                 val_axis_title = NULL,
                                 res_axis_title = NULL,
                                 FD_axis_title = NULL,
                                 aspectratio = c(2.5, .75, 1),
                                 showticklabels = TRUE) {
  if (!requireNamespace("lattice", quietly = TRUE)) {
    stop(
      "Package \"lattice\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # data <- X$data
  # curve_points <- X$fparam
  valRange <- c(
    floor(min(X$data,na.rm = TRUE)),
    ceiling(max(X$data,na.rm = TRUE))
  )

  name <- V1 <- value <- NULL
  data1 <- X$data
  colnames(data1) <- 1:ncol(X)
  plotData <- cbind(X$fparam,data1) %>%
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
  r_z <- range(plotData$Value,na.rm = TRUE)
  r_z[2] <- r_z[2]+4
  z_range <- round(r_z, -nchar(round(max(plotData$Value,na.rm = TRUE)))+1)

  if (showticklabels) {
    scale_info <- list(
      col = "black", arrows = FALSE, cex = 1.2,
      x = list(
        at = seq(min(X$fparam),
          max(X$fparam),
          length.out = 5
        ),
        # labels=.specify_decimal(seq(max(X$fparam),
        #                                 min(X$fparam),
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
        labels = c('',.specify_decimal(
          seq(z_range[1], z_range[2], length.out = 5)[-1], 0
        ))
      )
    )
  } else {
    scale_info <- list(
      col = "black", arrows = FALSE, cex = 0.75,
      x = list(
        at = seq(min(X$fparam),
          max(X$fparam),
          length.out = 5
        ),
        labels = NULL
      ),
      y = list(
        at = -seq(1, ncol(X), length.out = 7),
        labels = NULL
      ),
      z = list(
        at = seq(z_range[1], z_range[2],
          length.out = 5
        ),
        labels = NULL
      )
    )
  }

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
    xlim = rev(range(X$fparam)),
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
