#' Change Plot
#'
#' This plots the functional data with means for each segment.
#'
#' @inheritParams .plot_fd
#' @param warnings Boolean if warning should be given
#'
#' @return An interactive plotly plot for the data.
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_change(dfts(electricity$data[,1:150]),changes = c(66, 144))
#'    \item .plot_change(electricity,changes = c(66, 144, 204, 243, 305))
#'  }
.plot_change <- function(X, changes,
                  plot_title = NULL,
                  val_axis_title = "Value",
                  res_axis_title = "resolution",
                  FD_axis_title = "FD Sims",
                  eye = list(x = -1.5, y = -1.5, z = 1.5),
                  aspectratio = list(x = 1, y = 1, z = 1),
                  showticklabels = TRUE,
                  warnings = TRUE) {
  # Setup
  if(is.null(eye)) eye <- list(x = -1.5, y = -1.5, z = 1.5)
  if(is.null(aspectratio)) aspectratio <- list(x = 1, y = 1, z = 1)

  X <- dfts(X)
  data <- X$data
  curve_points <- X$fparam



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
  means <- rowMeans(data[, 1:min(changes),drop=FALSE])
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
  means <- rowMeans(data[, (max(changes) + 1):number,drop=FALSE])
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

    plotData_mean <- rbind(
      plotData_mean,
      data.frame(
        "resolution" = curve_points,
        "FDRep" = j,
        "Color" = length(changes) + 1,
        "Value" = means
      )
    )
  }

  # Color Additional Groups
  if (length(changes) > 1) {
    for (i in 2:length(changes)) {
      means <- rowMeans(data[, (changes[i - 1] + 1):changes[i], drop=FALSE])

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
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
  if (length(changes) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
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
    plotly::add_lines(data=plotData_mean,
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
}



#' Interval Plot
#'
#' A plot to show changes and their confidence interval on a functional data set.
#'
#' @inheritParams .plot_fd
#' @param intervals Data.frame from \code{confidence_interval()}.
#' @param highlight_changes Boolean if changes should be highlighted with a
#'  black line.
#'
#' @return An interactive plotly plot for the data.
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_interval(dfts(electricity$data[,1:80]),
#'                         data.frame('change'=66, 'lower'=60.17, 'upper'=72.54))
#'  }
.plot_interval <- function(X, intervals,
                          plot_title = X$name,
                          val_axis_title = "Value",
                          res_axis_title = "resolution",
                          FD_axis_title = "FD Sims",
                          eye = list(x = -1.5, y = -1.5, z = 1.5),
                          aspectratio = list(x = 1, y = 1, z = 1),
                          showticklabels = TRUE,
                          highlight_changes = TRUE, int.gradual=FALSE) {
  # Setup
  if(is.null(intervals)) stop('The variable `intervals` cannot be NULL.',
                              call. = FALSE)
  if(is.null(eye)) eye <- list(x = -1.5, y = -1.5, z = 1.5)
  if(is.null(aspectratio)) aspectratio <- list(x = 1, y = 1, z = 1)

  X <- dfts(X)
  if(!is.null(intervals)){
    changes <- unique(c(0,intervals[,1],ncol(X)))
    changes <- changes[order(changes)]

    # Get Colors
    plot_colors <- RColorBrewer::brewer.pal(
      min(9, max(length(changes) - 1,3)), "Set1")

    if (length(changes) > 9) {
      plot_colors <- rep(plot_colors,
                         ceiling(c(length(changes) + 1) / 9))[1:(length(changes) - 1)]
    }

  }else{
    changes <- 0:ncol(X)

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

  # Color and Group (Always at least 2 changes - start/end)
  for(i in 1:(length(changes)-1)){
    region <- (changes[i]+1):changes[i+1]
    means <- rowMeans(X$data[, region, drop=FALSE])
    for(j in region){
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$fparam,
          "FDRep" = j,#X$labels[j],
          "Color" = plot_colors[i],
          "Value" = X$data[, j]
        )
      )

      plotData_mean <- rbind(
        plotData_mean,
        data.frame(
          "resolution" = X$fparam,
          "FDRep" = j,#X$labels[j],
          "Color" = plot_colors[i],
          "Value" = means
        )
      )
    }

  }

  # Color Confidence Intervals and Changes
  if(!is.null(intervals)){

    if(!int.gradual){
      for(i in 1:nrow(intervals)){
        # Lower
        int <- floor(intervals[i,2]):(round(intervals[i,1]))
        plotData_mean[plotData_mean$FDRep %in% int, "Color"] <- 'gray'
        #X$labels[int], "Color"] <- 'gray'
        #colorRampPalette(c("lightgray", "darkgray"))(length(int))

        # Upper
        int <- (round(intervals[i,1])):round(intervals[i,3])
        plotData_mean[plotData_mean$FDRep %in% int, "Color"] <- 'gray'
                        #X$labels[int], "Color"] <- 'gray'

        # Changes
        if(highlight_changes){
          plotData_mean[plotData_mean$FDRep %in%
                          c(round(intervals$change[i]) + 0:1),
                        # X$labels[round(intervals$change[i]) + 0:1],
                        "Color"] <- 'black'

        }
      }

      if(highlight_changes){
        plot_colors <- c(plot_colors,'black','gray')
      }else{
        plot_colors <- c(plot_colors,'gray')
      }


      return_plot <- plotly::plot_ly(plotData,
                                     x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                                     type = "scatter3d", mode = "lines",
                                     split = ~ as.factor(FDRep),
                                     color = ~ as.factor(Color),
                                     colors = plot_colors,
                                     opacity=0.15
      ) %>%
        plotly::add_lines(data=plotData_mean,
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

    }else{

      for(i in 1:nrow(intervals)){
        change <- round(intervals$change[i])
        lower <- floor(intervals$lower[i])
        upper <- floor(intervals$upper[i])

        # Lower Ramp
        # rgb <- pmin(255,grDevices::col2rgb(plot_colors[i])*1.25)/255
        # light_plot_col <- grDevices::rgb(rgb[1],rgb[2],rgb[3])
        color_ramp <-
          grDevices::colorRampPalette(#c('darkgray','darkgray'))(
            c('darkgray',plot_colors[i]))(
            change - lower+2)
        if( change >= lower ){
          for (j in change:lower) {
            # Ensure not another break/interval
            curr_color <- plotData[plotData$FDRep==j,"Color"][1]
                                     #X$labels[j],"Color"][1]
            use_color <- ifelse(curr_color == plot_colors[i],
                                color_ramp[change-j+1],
                                min(curr_color, color_ramp[change-j+1]))

            # Save colors
            plotData[plotData$FDRep==j, "Color"] <- use_color
            plotData_mean[plotData_mean$FDRep==j, "Color"] <- use_color
            # plotData[plotData$FDRep==X$labels[j], "Color"] <- use_color
            # plotData_mean[plotData_mean$FDRep==X$labels[j], "Color"] <- use_color
          }
        }

        # # Upper Ramp
        color_ramp <-
          grDevices::colorRampPalette(
            c('darkgray',plot_colors[i+1]))(
              upper - change+2)
        if( (change+1) <= upper ){
          for (j in (change+1):upper) {
            # Ensure not another break/interval
            curr_color <- plotData[plotData$FDRep==j,"Color"][1]
                                     #X$labels[j],"Color"][1]
            use_color <- ifelse(curr_color == plot_colors[i+1],
                                color_ramp[j-(change+1)+1],
                                min(curr_color, color_ramp[j-(change+1)+1]))

            # Save colors
            plotData[plotData$FDRep==j, "Color"] <- use_color
            plotData_mean[plotData_mean$FDRep==j, "Color"] <- use_color
            # plotData[plotData$FDRep==X$labels[j], "Color"] <- use_color
            # plotData_mean[plotData_mean$FDRep==X$labels[j], "Color"] <- use_color
          }
        }

        if(highlight_changes){
          plotData_mean[plotData_mean$FDRep %in% c(change+0:1),"Color"] <- 'black'
          plotData[plotData$FDRep %in% c(change+0:1),"Color"] <- 'black'
          # plotData_mean[plotData_mean$FDRep %in% X$labels[change+0:1],"Color"] <- 'black'
          # plotData[plotData$FDRep %in% X$labels[change+0:1],"Color"] <- 'black'
        }
      }


      # Plot

      pal <- plotData_mean$Color
      names(pal) <- as.factor(plotData_mean$FDRep)
      # pal <- setNames(pal, as.factor(plotData_mean$FDRep))
      pal <- pal[unique(names(pal))]

      return_plot <- plotly::plot_ly(plotData,
                      x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                      type = "scatter3d", mode = "lines",
                      split = ~ as.factor(FDRep),
                      color = ~ as.factor(FDRep),
                      colors = pal,
                      opacity=0.15
      ) %>%
        plotly::add_lines(data=plotData_mean,
                          x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                          type = "scatter3d", mode = "lines",
                          split = ~ as.factor(FDRep),
                          color = ~ as.factor(FDRep),
                          colors = pal,#c(plot_colors,'black',ramp),
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
    }
  }

  return_plot
}
