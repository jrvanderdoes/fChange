#' Title
#'
#' @param data
#' @param curve_points
#' @param CPs
#' @param plot_title
#' @param val_axis_title
#' @param res_axis_title
#' @param FD_axis_title
#' @param eye
#' @param aspectratio
#' @param showticklabels
#'
#' @return
#' @export
#'
#' @examples
#' change_plot(electricity,curve_points = 1:nrow(electricity),CPs = c(66, 144, 204, 243, 305))
change_plot <- function(data, curve_points, CPs,
                  plot_title = NULL,
                  val_axis_title = "Value",
                  res_axis_title = "resolution",
                  FD_axis_title = "FD Sims",
                  eye = list(x = -1.5, y = -1.5, z = 1.5),
                  aspectratio = list(x = 1, y = 1, z = 1),
                  showticklabels = TRUE,
                  warnings = TRUE) {
  if(ncol(data)>750 & warnings)
    warning(paste0('Time to generate plot may be extensive ',
                   '(Ignore warning by setting warnings=FALSE).'),
            call. = FALSE)

  number <- length(data[1, ])
  valRange <- c(min(data), max(data))

  plotData <- plotData_mean <- data.frame(
    "resolution" = NA,
    "FDRep" = NA,
    "Color" = NA,
    "Value" = NA
  )[-1, ]

  # Color Group to first CP
  means <- rowMeans(data[, 1:min(CPs),drop=FALSE])
  for (j in 1:min(CPs)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = 1,
        "Value" = data[, j]
      )
    )

    plotData_mean <- rbind(
      plotData_mean,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = 1,
        "Value" = means
      )
    )
  }
  # Color Group from last CP
  means <- rowMeans(data[, (max(CPs) + 1):number,drop=FALSE])
  for (j in (max(CPs) + 1):number) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = length(CPs) + 1,
        "Value" = data[, j]
      )
    )

    plotData_mean <- rbind(
      plotData_mean,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = length(CPs) + 1,
        "Value" = means
      )
    )
  }

  # Color Additional Groups
  if (length(CPs) > 1) {
    for (i in 2:length(CPs)) {
      means <- rowMeans(data[, (CPs[i - 1] + 1):CPs[i], drop=FALSE])

      for (j in (CPs[i - 1] + 1):CPs[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = curve_points,
            "FDRep" = j,
            "Color" = i,
            "Value" = data[, j]
          )
        )

        plotData_mean <- rbind(
          plotData_mean,
          data.frame(
            "resolution" = curve_points,
            "FDRep" = j,
            "Color" = i,
            "Value" = means
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

  plotly::plot_ly(plotData,
                  x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  split = ~ as.factor(FDRep),
                  # color = 'gray',
                  # colors='gray',
                  color = ~ as.factor(Color),
                  colors = tmpColors,
                  opacity=0.15
  ) %>%
    plotly::add_lines(., data=plotData_mean,
                    x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                    type = "scatter3d", mode = "lines",
                    split = ~ as.factor(FDRep),
                    color = ~ as.factor(Color),
                    colors = tmpColors,
                    opacity=1
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


#' #' Title
#' #'
#' #' @param X
#' #' @param intervals
#' #' @param plot_title
#' #' @param val_axis_title
#' #' @param res_axis_title
#' #' @param FD_axis_title
#' #' @param eye
#' #' @param aspectratio
#' #' @param showticklabels
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' interval_plot <- function(X, intervals,
#'                         plot_title = X$name,
#'                         val_axis_title = "Value",
#'                         res_axis_title = "resolution",
#'                         FD_axis_title = "FD Sims",
#'                         eye = list(x = -1.5, y = -1.5, z = 1.5),
#'                         aspectratio = list(x = 1, y = 1, z = 1),
#'                         showticklabels = TRUE) {
#'   X <- .check_data(X)
#'
#'   data <- X$data
#'   curve_points <- X$intraobs
#'   CPs <- intervals[,1]
#'
#'   number <- length(data[1, ])
#'   valRange <- c(min(data), max(data))
#'
#'
#'   plotData <- plotData_mean <- data.frame(
#'     "resolution" = NA,
#'     "FDRep" = NA,
#'     "Color" = NA,
#'     "Value" = NA
#'   )[-1, ]
#'
#'   # Color Group to first CP
#'   means <- rowMeans(data[, 1:min(CPs),drop=FALSE])
#'   for (j in 1:min(CPs)) {
#'     plotData <- rbind(
#'       plotData,
#'       data.frame(
#'         "resolution" = curve_points,
#'         "FDRep" = j,
#'         "Color" = 1,
#'         "Value" = data[, j]
#'       )
#'     )
#'
#'     plotData_mean <- rbind(
#'       plotData_mean,
#'       data.frame(
#'         "resolution" = curve_points,
#'         "FDRep" = j,
#'         "Color" = 1,
#'         "Value" = means
#'       )
#'     )
#'   }
#'   # Color Group from last CP
#'   means <- rowMeans(data[, (max(CPs) + 1):number,drop=FALSE])
#'   for (j in (max(CPs) + 1):number) {
#'     plotData <- rbind(
#'       plotData,
#'       data.frame(
#'         "resolution" = curve_points,
#'         "FDRep" = j,
#'         "Color" = length(CPs) + 1,
#'         "Value" = data[, j]
#'       )
#'     )
#'
#'     plotData_mean <- rbind(
#'       plotData_mean,
#'       data.frame(
#'         "resolution" = curve_points,
#'         "FDRep" = j,
#'         "Color" = length(CPs) + 1,
#'         "Value" = means
#'       )
#'     )
#'   }
#'
#'   # Color Additional Groups
#'   if (length(CPs) > 1) {
#'     for (i in 2:length(CPs)) {
#'       means <- rowMeans(data[, (CPs[i - 1] + 1):CPs[i], drop=FALSE])
#'
#'       for (j in (CPs[i - 1] + 1):CPs[i]) {
#'         plotData <- rbind(
#'           plotData,
#'           data.frame(
#'             "resolution" = curve_points,
#'             "FDRep" = j,
#'             "Color" = i,
#'             "Value" = data[, j]
#'           )
#'         )
#'
#'         plotData_mean <- rbind(
#'           plotData_mean,
#'           data.frame(
#'             "resolution" = curve_points,
#'             "FDRep" = j,
#'             "Color" = i,
#'             "Value" = means
#'           )
#'         )
#'       }
#'     }
#'   }
#'
#'   # Color Confidence Intervals
#'   for(i in 1:nrow(intervals)){
#'     # Lower
#'     plotData_mean[plotData_mean$FDRep %in%
#'                     floor(intervals[i,2]):round(intervals[i,1]),"Color"] <- 'gray'
#'
#'     # Upper
#'     plotData_mean[plotData_mean$FDRep %in%
#'                     round(intervals[i,1]):floor(intervals[i,3]),"Color"] <- 'gray'
#'   }
#'
#'   scene <- list(
#'     camera = list(eye = eye),
#'     aspectmode = "manual",
#'     aspectratio = aspectratio
#'   )
#'
#'   # Get Colors
#'   tmpColors <- RColorBrewer::brewer.pal(length(CPs) + 1, "Set1")
#'   tmpColors <- c(tmpColors,'gray')
#'
#'   # Parameters
#'   if(is.null(res_axis_title)){
#'     res_axis_title <- ''
#'   }
#'   if(is.null(FD_axis_title)){
#'     FD_axis_title <- ''
#'   }
#'   if(is.null(val_axis_title)){
#'     val_axis_title <- ''
#'   }
#'
#'
#'   plotly::plot_ly(plotData,
#'                   x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
#'                   type = "scatter3d", mode = "lines",
#'                   split = ~ as.factor(FDRep),
#'                   # color = 'gray',
#'                   # colors='gray',
#'                   color = ~ as.factor(Color),
#'                   colors = tmpColors,
#'                   opacity=0.15
#'   ) %>%
#'     plotly::add_lines(., data=plotData_mean,
#'                       x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
#'                       type = "scatter3d", mode = "lines",
#'                       split = ~ as.factor(FDRep),
#'                       color = ~ as.factor(Color),
#'                       colors = tmpColors,
#'                       opacity=1
#'     ) %>%
#'     plotly::layout(
#'       scene = list(
#'         yaxis = list(
#'           title = res_axis_title,
#'           showticklabels = showticklabels
#'         ),
#'         xaxis = list(
#'           title = FD_axis_title,
#'           showticklabels = showticklabels
#'         ),
#'         zaxis = list(
#'           title = val_axis_title,
#'           showticklabels = showticklabels
#'         )
#'       )
#'     ) %>%
#'     plotly::layout(title = plot_title, scene = scene) %>%
#'     plotly::layout(showlegend = FALSE)
#' }

#' Title
#'
#' @param X
#' @param intervals
#' @param plot_title
#' @param val_axis_title
#' @param res_axis_title
#' @param FD_axis_title
#' @param eye
#' @param aspectratio
#' @param showticklabels
#'
#' @return
#' @export
#'
#' @examples
interval_plot <- function(X, intervals,
                          plot_title = X$name,
                          val_axis_title = "Value",
                          res_axis_title = "resolution",
                          FD_axis_title = "FD Sims",
                          eye = list(x = -1.5, y = -1.5, z = 1.5),
                          aspectratio = list(x = 1, y = 1, z = 1),
                          showticklabels = TRUE,
                          highlight_changes = TRUE) {
  X <- .check_data(X)
  if(!is.null(intervals)){
    CPs <- unique(c(0,intervals[,1],ncol(X)))
    CPs <- CPs[order(CPs)]

    # Get Colors
    plot_colors <- RColorBrewer::brewer.pal(length(CPs) - 1, "Set1")
  }else{
    CPs <- 0:ncol(X)

    plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
    plot_colors[6] <- "yellow"
  }


  # Plot Parameters
  plotData <- plotData_mean <- data.frame()

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  if(is.null(res_axis_title)){
    res_axis_title <- ''
  }
  if(is.null(FD_axis_title)){
    FD_axis_title <- ''
  }
  if(is.null(val_axis_title)){
    val_axis_title <- ''
  }

  # Color and Group (Always at least 2 CPs - start/end)
  for(i in 1:(length(CPs)-1)){
    region <- (CPs[i]+1):CPs[i+1]
    means <- rowMeans(X$data[, region, drop=FALSE])
    for(j in region){
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$intraobs,
          "FDRep" = X$labels[j],
          "Color" = plot_colors[i],
          "Value" = X$data[, j]
        )
      )

      plotData_mean <- rbind(
        plotData_mean,
        data.frame(
          "resolution" = X$intraobs,
          "FDRep" = X$labels[j],
          "Color" = plot_colors[i],
          "Value" = means
        )
      )
    }

  }

  # Color Confidence Intervals and Changes
  if(!is.null(intervals)){

    for(i in 1:nrow(intervals)){
      # Lower
      int <- floor(intervals[i,2]):(round(intervals[i,1]))
      plotData_mean[plotData_mean$FDRep %in% X$labels[int], "Color"] <- 'gray'
      #colorRampPalette(c("lightgray", "darkgray"))(length(int))

      # Upper
      int <- (round(intervals[i,1])):round(intervals[i,3])
      plotData_mean[plotData_mean$FDRep %in% X$labels[int], "Color"] <- 'gray'

      # Changes
      if(highlight_changes){
        plotData_mean[plotData_mean$FDRep %in%
                        X$labels[round(intervals$change[i]) + 0:1],
                      "Color"] <- 'black'

      }
    }

    if(highlight_changes){
      plot_colors <- c(plot_colors,'black','gray')
    }else{
      plot_colors <- c(plot_colors,'gray')
    }
  }

  # Plot
  plotly::plot_ly(plotData,
                  x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  split = ~ as.factor(FDRep),
                  color = ~ as.factor(Color),
                  colors = plot_colors,
                  opacity=0.15
  ) %>%
    plotly::add_lines(., data=plotData_mean,
                      x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                      type = "scatter3d", mode = "lines",
                      split = ~ as.factor(FDRep),
                      color = ~ as.factor(Color),
                      colors = plot_colors,
                      opacity=1
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
