
#' Plot functional data
#'
#' \code{plot_fd} plots functional data either in an fd object or evaluated at certain points
#'
#' @param data An fd object or a data.frame of evaluated fd objects
#'     (the columns bieng lines and the rows the evaluated points)
#'
#' @param eval_points An optional vector containing the points at which the
#'     evaluated data was evaluated at
#' @param plot_title Optional string to title the plot
#'
#' @return A plotly plot
#'
#' @export
#'
#' @examples
#' # Generate Data with a mean change point
#' evalPts <- seq(0,1,0.05)
#' data_KL <- generate_data_fd(ns = c(100,100),
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
#' # Plot an FD object. Note differences are due to plotting the smoothed data
#' data_fd <- fda::Data2fd(argvals=evalPts, y=as.matrix(data_KL),
#'                    basisobj=fda::create.bspline.basis(rangeval=range(evalPts)))
#' plot_fd(data=data_fd)
plot_fd <- function(data, eval_points=NULL, plot_title=NULL){

  if(fda::is.fd(data)){
    eval_Range <- data$basis$rangeval
    eval_points <- seq(eval_Range[1],eval_Range[2],length.out=100)

    fd_eval <- fda::eval.fd(eval_points, data)

    plot <- .plot_evalfd_3dlines_plotly(fd_eval, eval_points, plot_title)
  }else if(!is.null(eval_points)){
    plot <- .plot_evalfd_3dlines_plotly(data, eval_points, plot_title)
  }else{
    stop('Either data must be an fd object or eval_points must exist for the evaluated data')
  }

  plot
}


.plot_evalfd_3dlines_plotly <- function(fd_eval, singleRangeSeq, titleVal=NULL){

  number <- length(fd_eval[1,])
  valRange <- c(min(fd_eval),max(fd_eval))

  plotData <- data.frame('evalRange'=singleRangeSeq,
                         'FDRep'=1,
                         'Value'=fd_eval[,1])

  for(i in 2:number){
    plotData <- rbind(plotData,
                      data.frame('evalRange'=singleRangeSeq,
                                 'FDRep'=i,
                                 'Value'=fd_eval[,i])
    )
  }

  scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))

  tmpColors <- RColorBrewer::brewer.pal(11,"Spectral")
  tmpColors[6] <- 'yellow'

  magrittr::`%>%`(
    magrittr::`%>%`(
      magrittr::`%>%`(
        plotly::plot_ly(plotData,
                        x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
                        type = 'scatter3d', mode = 'lines',
                        color = ~as.factor(FDRep),
                        colors = tmpColors),
        plotly::layout(
          scene = list(
            yaxis = list(title = "EvalRange"),
            xaxis = list(title = "FD Sims"),
            zaxis = list(title = "Value")
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
