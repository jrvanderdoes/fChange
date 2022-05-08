
#' Plot functional data
#'
#' \code{plot_fd} plots functional data either in an fd object or evaluated at certain points
#'
#' @param data An fd object or a data.frame of evaluated fd objects
#'     (the columns being lines and the rows the evaluated points)
#' @param eval_points An vector containing the points at which the
#'     evaluated data was evaluated at
#' @param CPs (Optional) Vectors of numeric values indicating the location of
#'     change points
#' @param plot_title Optional string to title the plot
#'
#' @return A plotly plot
#'
#' @export
#'
#' @examples
#' # Generate Data with a mean change point
#' evalPts <- seq(0,1,0.05)
#' data_KL <- generate_data_fd(ns = c(200,200),
#'                                eigsList = list(c(3,2,1,0.5),
#'                                                c(3,2,1,0.5)),
#'                                basesList = list(fda::create.bspline.basis(nbasis=4, norder=4),
#'                                                 fda::create.bspline.basis(nbasis=4, norder=4)),
#'                                meansList = c(-1,1),
#'                                distsArray = c('Normal','Normal'),
#'                                evals = evalPts,
#'                                kappasArray = c(0,0))
#' plot_fd(data=data_KL,eval_points=evalPts)
#'
#' plot_fd(data=data_KL,eval_points=evalPts, CPs=200)
#' plot_fd(data=data_KL,eval_points=evalPts, CPs=seq(40,360,20))
plot_fd <- function(data, eval_points = 1:nrow(data), CPs=NULL, plot_title=NULL,
                    val_axis_title = 'Value', eval_axis_title='EvalRange',
                    FD_axis_title = 'FD Sims',FDReps=1:ncol(data)){

  if(!is.null(CPs)){
    plot <- .plot_evalfd_3dlines_cps(data, eval_points, CPs, plot_title,
                                     val_axis_title, eval_axis_title,
                                     FD_axis_title)#,FDReps)
  }else{
    plot <- .plot_evalfd_3dlines(data, eval_points, plot_title,
                                 val_axis_title, eval_axis_title,
                                 FD_axis_title,FDReps)
  }

  plot
}


.plot_evalfd_3dlines <- function(fd_eval, singleRangeSeq, titleVal=NULL,
                                 val_axis_title = 'Value',
                                 eval_axis_title='EvalRange',
                                 FD_axis_title = 'FD Sims',
                                 FDReps = 1:ncol(fd_eval)){

  number <- length(fd_eval[1,])
  valRange <- c(min(fd_eval),max(fd_eval))

  plotData <- data.frame('evalRange'=singleRangeSeq,
                         'FDRep'=FDReps[1],
                         'Value'=fd_eval[,1])

  for(i in 2:number){
    plotData <- rbind(plotData,
                      data.frame('evalRange'=singleRangeSeq,
                                 'FDRep'=FDReps[i],
                                 'Value'=as.factor(fd_eval[,i]))
    )
  }

  scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))

  tmpColors <- RColorBrewer::brewer.pal(11,"Spectral")
  tmpColors[6] <- 'yellow'

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~FDRep, y = ~evalRange, z = ~Value,
                        type = 'scatter3d', mode = 'lines',
                        color = ~ as.factor(FDRep),
                        colors = tmpColors),
        plotly::layout(
          scene = list(
            yaxis = list(title = eval_axis_title),
            xaxis = list(title = FD_axis_title),
            zaxis = list(title = val_axis_title)
          ))),
      plotly::layout(title = titleVal, scene = scene)),
    plotly::layout(showlegend = FALSE))

  #plotly::plot_ly(plotData,
  #                x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
  #                type = 'scatter3d', mode = 'lines',
  #                color = ~as.factor(FDRep),
  #                colors = tmpColors) %>%
  #  #colors = c("red", "yellow", "blue")) %>%
  #  #colors='Spectral') %>%
  #  plotly::layout(
  #    scene = list(
  #      yaxis = list(title = "EvalRange"),
  #      xaxis = list(title = "FD Sims"),
  #      zaxis = list(title = "Value")
  #    )) %>%
  #  plotly::layout(title = titleVal, scene = scene) %>%
  #  plotly::layout(showlegend = FALSE)
  ##line = list(width = 4, color = ~as.factor(FDRep),
  ##            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
}

.plot_evalfd_3dlines_cps <- function(fd_eval, singleRangeSeq, CPs,
                                     titleVal=NULL,
                                     val_axis_title = 'Value',
                                     eval_axis_title='EvalRange',
                                     FD_axis_title = 'FD Sims'){

  number <- length(fd_eval[1,])
  valRange <- c(min(fd_eval),max(fd_eval))

  plotData <- data.frame('evalRange'=NA,
                         'FDRep'=NA,
                         'Color'=NA,
                         'Value'=NA)[-1,]

  # Color Group to first CP
  for(j in 1:min(CPs)){
    plotData <- rbind(plotData,
                      data.frame('evalRange'=singleRangeSeq,
                                 'FDRep'=j,
                                 'Color'=1,
                                 'Value'=fd_eval[,j]))
  }
  # Color Group from last CP
  for(j in (max(CPs)+1):number){
    plotData <- rbind(plotData,
                      data.frame('evalRange'=singleRangeSeq,
                                 'FDRep'=j,
                                 'Color'=length(CPs)+1,
                                 'Value'=fd_eval[,j]))
  }

  # Color Additional Groups
  if(length(CPs)>1){
    for(i in 2:length(CPs)){
      for(j in (CPs[i-1]+1):CPs[i])
        plotData <- rbind(plotData,
                          data.frame('evalRange'=singleRangeSeq,
                                     'FDRep'=j,
                                     'Color'=i,
                                     'Value'=fd_eval[,j]))
    }
  }

  scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9,max(3,length(CPs)+1)),"Set1")
  if(length(CPs)>9)
    tmpColors <- rep(tmpColors, ceiling(c(length(CPs)+1)/9))[1:(length(CPs)+1)]

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
                        type = 'scatter3d', mode = 'lines',
                        split=~as.factor(FDRep),
                        color = ~as.factor(Color),
                        colors = tmpColors),
        plotly::layout(
          scene = list(
            yaxis = list(title = eval_axis_title),
            xaxis = list(title = FD_axis_title),
            zaxis = list(title = val_axis_title)
          ))),
      plotly::layout(title = titleVal, scene = scene)),
    plotly::layout(showlegend = FALSE))
}

