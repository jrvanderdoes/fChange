
#' Plot functional data
#'
#' \code{plot_fd} plots functional data either in an fd object or evaluated at certain points
#'
#' @param data A data.frame of evaluated fd objects (the columns being lines and
#'     the rows the evaluated points)
#' @param CPs (Optional) Vectors of numeric values indicating the location of
#'     change points. This will color each section differently. Default vallue
#'     is NULL.
#' @param curve_points (Optional) An vector containing the points at which the
#'     data was evaluated. Default is 1:nrow(data)
#' @param plot_title (Optional) String to title the plot. Default value is NULL.
#' @param val_axis_title (Optional) String to title the axis for values of
#'     observation. Default value is 'Value'.
#' @param res_axis_title (Optional) String to title the axis for the resolution
#'     range axis (observations in an FD object). Default value is 'resolution'.
#' @param FD_axis_title (Optional) String to title the axis for FD observations.
#'     Default value is 'Observation'.
#' @param FDReps (Optional) Vector of values for FD observation names. Default
#'     value is 1:ncol(data).
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
#' @export
#'
#' @examples
#' plot_fd(data=electricity[,1:100], CPs=c(50))
#' plot_fd(data=electricity, CPs=c(50,150,220,300),
#'         interactive=FALSE, showticklabels=FALSE)
plot_fd <- function(data, CPs=NULL, curve_points = 1:nrow(data), plot_title=NULL,
                    val_axis_title = 'Value', res_axis_title='resolution',
                    FD_axis_title = 'Observations',FDReps=1:ncol(data),
                    eye = list(x = -1.5, y = -1.5, z = 1.5),
                    aspectratio=list(x=1,y=1,z=1),
                    showticklabels=TRUE, interactive=TRUE){

  if(!interactive){
    fdPlot <- .plot_evalfd_highdim(data=data, curve_points=curve_points,
                                  CPs=CPs, plot_title=plot_title,
                                  val_axis_title=val_axis_title,
                                  res_axis_title=res_axis_title,
                                  FD_axis_title=FD_axis_title,
                                  aspectratio=aspectratio,
                                  showticklabels=showticklabels)

  } else if(!is.null(CPs) && length(stats::na.omit(CPs))>0){
    fdPlot <- .plot_evalfd_3dlines_cps(data=data, curve_points=curve_points,
                                       CPs=CPs[order(CPs)], plot_title=plot_title,
                                       val_axis_title=val_axis_title,
                                       res_axis_title=res_axis_title,
                                       FD_axis_title=FD_axis_title,
                                       eye=eye, aspectratio=aspectratio,
                                       showticklabels=showticklabels)#,FDReps)
  }else{
    fdPlot <- .plot_evalfd_3dlines(data=data, curve_points=curve_points,
                                   plot_title=plot_title,
                                   val_axis_title=val_axis_title,
                                   res_axis_title=res_axis_title,
                                   FD_axis_title=FD_axis_title,FDReps=FDReps,
                                   eye=eye, aspectratio=aspectratio,
                                   showticklabels=showticklabels)
  }

  fdPlot
}


#' Plot Functional Data Without Change Points
#'
#' This (internal) function to plot the function data, with no coloring based on
#'     change points.
#'
#' @inheritParams plot_fd
#'
#' @return A plotly plot
#'
#' @noRd
.plot_evalfd_3dlines <- function(data, curve_points, plot_title=NULL,
                                 val_axis_title = 'Value',
                                 res_axis_title='Resolution',
                                 FD_axis_title = 'Observation',
                                 FDReps = 1:ncol(data),
                                 eye = list(x = -1.5, y = -1.5, z = 1.5),
                                 aspectratio=list(x=1,y=1,z=1),
                                 showticklabels=TRUE){

  number <- length(data[1,])
  valRange <- c(min(data),max(data))

  plotData <- data.frame('resolution'=curve_points,
                         'FDRep'=FDReps[1],
                         'Value'=data[,1])

  for(i in 2:number){
    plotData <- rbind(plotData,
                      data.frame('resolution'=curve_points,
                                 'FDRep'=FDReps[i],
                                 'Value'=as.factor(data[,i]))
    )
  }

  scene <- list(camera = list(eye = eye),
                             aspectmode = "manual",
                             aspectratio = aspectratio)

  tmpColors <- RColorBrewer::brewer.pal(11,"Spectral")
  tmpColors[6] <- 'yellow'

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~FDRep, y = ~resolution, z = ~Value,
                        type = 'scatter3d', mode = 'lines',
                        color = ~ as.factor(FDRep),
                        colors = tmpColors),
        plotly::layout(
          scene = list(
            yaxis = list(title = res_axis_title,
                         showticklabels=showticklabels),
            xaxis = list(title = FD_axis_title,
                         showticklabels=showticklabels),
            zaxis = list(title = val_axis_title,
                         showticklabels=showticklabels)
          ))),
      plotly::layout(title = plot_title, scene = scene)),
    plotly::layout(showlegend = FALSE))

  #plotly::plot_ly(plotData,
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
  ##line = list(width = 4, color = ~as.factor(FDRep),
  ##            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
}


#' Plot Functional Data With Change Points
#'
#' This (internal) function to plot the function data, with coloring based on
#'     change points.
#'
#' @inheritParams plot_fd
#' @param CPs Vectors of numeric values indicating the location of change points.
#'     This will color each section differently.
#'
#' @return A plotly plot
#'
#' @noRd
.plot_evalfd_3dlines_cps <- function(data, curve_points, CPs,
                                     plot_title=NULL,
                                     val_axis_title = 'Value',
                                     res_axis_title='resolution',
                                     FD_axis_title = 'FD Sims',
                                     eye = list(x = -1.5, y = -1.5, z = 1.5),
                                     aspectratio=list(x=1,y=1,z=1),
                                     showticklabels=TRUE){

  number <- length(data[1,])
  valRange <- c(min(data),max(data))

  plotData <- data.frame('resolution'=NA,
                         'FDRep'=NA,
                         'Color'=NA,
                         'Value'=NA)[-1,]

  # Color Group to first CP
  for(j in 1:min(CPs)){
    plotData <- rbind(plotData,
                      data.frame('resolution'=curve_points,
                                 'FDRep'=j,
                                 'Color'=1,
                                 'Value'=data[,j]))
  }
  # Color Group from last CP
  for(j in (max(CPs)+1):number){
    plotData <- rbind(plotData,
                      data.frame('resolution'=curve_points,
                                 'FDRep'=j,
                                 'Color'=length(CPs)+1,
                                 'Value'=data[,j]))
  }

  # Color Additional Groups
  if(length(CPs)>1){
    for(i in 2:length(CPs)){
      for(j in (CPs[i-1]+1):CPs[i])
        plotData <- rbind(plotData,
                          data.frame('resolution'=curve_points,
                                     'FDRep'=j,
                                     'Color'=i,
                                     'Value'=data[,j]))
    }
  }

  scene <- list(camera = list(eye = eye),
                aspectmode = "manual",
                aspectratio = aspectratio)

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9,max(3,length(CPs)+1)),"Set1")
  if(length(CPs)>9)
    tmpColors <- rep(tmpColors, ceiling(c(length(CPs)+1)/9))[1:(length(CPs)+1)]

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~as.factor(FDRep), y = ~resolution, z = ~Value,
                        type = 'scatter3d', mode = 'lines',
                        split=~as.factor(FDRep),
                        color = ~as.factor(Color),
                        colors = tmpColors),
        plotly::layout(
          scene = list(
            yaxis = list(title = res_axis_title,
                         showticklabels=showticklabels),
            xaxis = list(title = FD_axis_title,
                         showticklabels=showticklabels),
            zaxis = list(title = val_axis_title,
                         showticklabels=showticklabels)
          ))),
      plotly::layout(title = plot_title, scene = scene)),
    plotly::layout(showlegend = FALSE))
}


#' Plot With Surface for Speed
#'
#' @inheritParams plot_fd
#' @param aspectratio (Optional) List with certain parameters to determine the
#'     image size parameters. The default value is c(2.5,.75,1). Also any lists
#'     are converted to this (as it is likely from default call `in plot_fd()`).
#'
#' @return A lattice plot
#'
#' @noRd
.plot_evalfd_highdim <- function(data, curve_points, CPs=NULL,
                                plot_title=NULL,val_axis_title=NULL,
                                res_axis_title=NULL,
                                FD_axis_title=NULL,
                                aspectratio=c(2.5,.75,1),
                                showticklabels=FALSE){
  if (!requireNamespace("lattice", quietly = TRUE)) {
    stop(
      "Package \"lattice\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # This resets ratio, since default was likely used
  if(is.list(aspectratio)) aspectratio <- c(2.5,0.75,1)

  number <- length(data[1,])
  valRange <- c(floor(min(data)),
                ceiling(max(data)))

  plotData <- data.frame('resolution'=curve_points,
                         'FDRep'=1,
                         'Value'=data[,1])

  for(i in 2:number){
    plotData <- rbind(plotData,
                      data.frame('resolution'=curve_points,
                                 'FDRep'=i,
                                 'Value'=data[,i])
    )
  }

  ## Setup up Colors
  #   Rainbow for no CPs, colored for CPs
  plotData[['color']] <- rep(1:ncol(data), each=nrow(data))
  if(!is.null(CPs)){
    tmp_colors <- RColorBrewer::brewer.pal(min(9,max(3,length(CPs)+1)),"Set1")
    if(length(CPs)>9)
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs)+1)/9))[1:(length(CPs)+1)]

    CPs <- unique(c(1,CPs,ncol(data)))
    colors_plot <- rep(tmp_colors[1], ncol(data))
    for(i in 2:(length(CPs)-1)){
      colors_plot[CPs[i]:CPs[i+1]] <- tmp_colors[i]
    }
  } else{
    colors_plot <- RColorBrewer::brewer.pal(11,"Spectral")
    colors_plot[6] <- 'yellow'
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(data))
    plotData[['color']] <- rep(colors_plot,
                               each=nrow(data))
  }

  ## Set up Tick Labels
  z_range <- round(range(plotData$Value),-2)
  if(showticklabels){
    scale_info <- list(
      col='black',arrows=FALSE, cex= 0.75, cex.title=1.5,
      x = list(
        at=seq(min(curve_points),
               max(curve_points),
               length.out=5),
        # labels=.specify_decimal(rev(seq(max(curve_points),
        #                                 min(curve_points),
        #                                 length.out=4)),2)),
        labels=NULL),
      y = list(
        at=-seq(1, number, length.out=8),
        labels=.specify_decimal(
          seq(1, number, length.out=8),0)),
      z = list(
        at=seq(z_range[1], z_range[2],
                length.out=5),
        labels=.specify_decimal(
          seq(z_range[1], z_range[2], length.out=5),0)) )
  } else{
    scale_info <- list(
      col='black',arrows=FALSE, cex= 0.75, cex.title=1.5,
      x = list(
        at=seq(min(curve_points),
               max(curve_points),
               length.out=5),
        labels=NULL),
      y = list(
        at=-seq(1, number, length.out=8),
        labels=NULL),
      z = list(
        at=seq(z_range[1], z_range[2],
               length.out=5),
        labels=NULL) )
  }

  color <- NULL # Check note removal
  #lattice:::wireframe
  lattice::cloud(x=Value ~ (resolution)*(-FDRep),
                     data= plotData,
                     type="l",groups=color,
                     par.box = c(col = "transparent"),
                     par.settings =
                      list(axis.line=list(col="transparent"),
                           superpose.line = list(col=colors_plot)),
                     #screen=list(z = 90, x = -75,y=-45),
                     #trellis.par.set(list(axis.text=list(cex=2)),
                     #                "axis.line",list(col=NA)),
                     zlim=valRange,
                     aspect=aspectratio,
                     drape=TRUE,colorkey = FALSE,
                     scales = scale_info,
                     xlab=list(res_axis_title, rot = 30),
                     ylab=list(paste0("\n",FD_axis_title), rot = -30),
                     zlab=list(val_axis_title, rot = 90,just=0.75),
                     main=plot_title)
}
