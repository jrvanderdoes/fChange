
#########################
# .plot_fd_3dsurf <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- eval.fd(singleRangeSeq,
#                      data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   #trellis.par.set("axis.line",list(col='black'))
#   lattice::trellis.par.set("axis.line",list(col=NA))
#   lattice::wireframe(x=Value ~ (evalRange)*(-FDRep),
#                      data= plotData,
#                      #trellis.par.set(list(axis.text=list(cex=2)),
#                      #                "axis.line",list(col=NA)),
#                      zlim=valRange,
#                      aspect=c(3,.75,1),
#                      drape=TRUE,colorkey = FALSE,
#                      scales = list(arrows=FALSE,cex= 0.75,
#                                    cex.title=1.5,
#                                    x = list(at=seq(singleRange[1],
#                                                    singleRange[2],
#                                                    0.2),
#                                             labels=rev(seq(singleRange[2],
#                                                            singleRange[1],
#                                                            -0.2))),
#                                    y = list(at=-seq(1, number, 9),
#                                             labels=seq(1, number, 9))),
#                      xlab="Eval Range",ylab="\nFD Reps",zlab="Value",
#                      main=titleVal)
# }
#
#########################
# .plot_fd_3dsurf_plotly <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- eval.fd(singleRangeSeq,
#                      data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   plotly::plot_ly() %>%
#     plotly::add_trace(data = plotData,
#                       x=plotData$FDRep,
#                       y=plotData$evalRange,
#                       z=plotData$Value, type="mesh3d" )
#
#   tmpData <-
#     as.matrix(plotData %>%
#                 tidyr::pivot_wider(., names_from ='evalRange',values_from= 'Value')
#     )[,-1]
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   plotly::plot_ly(x = 1:number,
#                   y = singleRangeSeq,
#                   z = t(tmpData)) %>%
#     plotly::add_surface() %>%
#     plotly::layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     plotly::layout(title = titleVal, scene = scene)
# }
#
#
#########################
# .plot_evalfd_3dsurf <- function(fd_eval, singleRangeSeq,
#                                 titleVal=NULL){
#
#   number <- length(fd_eval[1,])
#   valRange <- c(floor(min(fd_eval)),
#                 ceiling(max(fd_eval)))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   #trellis.par.set("axis.line",list(col='black'))
#   withr::with_locale(lattice::trellis.par.set("axis.line",list(col=NA)),
#                      plot <-  lattice::wireframe(x=Value ~ (evalRange)*(-FDRep),
#                                                  data= plotData,
#                                                  #trellis.par.set(list(axis.text=list(cex=2)),
#                                                  #                "axis.line",list(col=NA)),
#                                                  zlim=valRange,
#                                                  aspect=c(3,.75,1),
#                                                  drape=TRUE,colorkey = FALSE,
#                                                  scales = list(arrows=FALSE,cex= 0.75,
#                                                                cex.title=1.5,
#                                                                x = list(at=seq(min(singleRangeSeq),
#                                                                                max(singleRangeSeq),
#                                                                                0.2),
#                                                                         labels=rev(seq(max(singleRangeSeq),
#                                                                                        min(singleRangeSeq),
#                                                                                        -0.2))),
#                                                                y = list(at=-seq(1, number, 9),
#                                                                         labels=seq(1, number, 9))),
#                                                  xlab="Eval Range",ylab="\nFD Reps",zlab="Value",
#                                                  main=titleVal)
#   )
#
#   plot
# }
#
#########################
# .plot_evalfd_3dsurf_plotly <- function(fd_eval, singleRangeSeq,
#                                        titleVal=NULL){
#
#   number <- length(fd_eval[1,])
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   plotly::plot_ly() %>%
#     plotly::add_trace(data = plotData,
#                       x=plotData$FDRep,
#                       y=plotData$evalRange,
#                       z=plotData$Value, type="mesh3d" )
#
#   tmpData <-
#     as.matrix(plotData %>%
#                 tidyr::pivot_wider(., names_from ='evalRange',values_from= 'Value')
#     )[,-1]
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   plotly::plot_ly(x = 1:number,
#                   y = singleRangeSeq,
#                   z = t(tmpData)) %>%
#     plotly::add_surface() %>%
#     plotly::layout(
#       scene = list(
#         yaxis = list(title = "EvalRange"),
#         xaxis = list(title = "FD Sims"),
#         zaxis = list(title = "Value")
#       )) %>%
#     plotly::layout(title = titleVal, scene = scene)
# }
#
#########################
# .plot_fd_3dlines_plotly <- function(fdobj, titleVal){
#
#   singleRange <- fdobj$basis$rangeval
#   singleRangeSeq <- seq(singleRange[1],singleRange[2],length.out=100)
#   number <- length(fdobj$fdnames$reps)
#
#   fd_eval <- fda::eval.fd(singleRangeSeq,
#                           data_fd)
#   valRange <- c(min(fd_eval),max(fd_eval))
#
#   plotData <- data.frame('evalRange'=singleRangeSeq,
#                          'FDRep'=1,
#                          'Value'=fd_eval[,1])
#
#   for(i in 2:number){
#     plotData <- rbind(plotData,
#                       data.frame('evalRange'=singleRangeSeq,
#                                  'FDRep'=i,
#                                  'Value'=fd_eval[,i])
#     )
#   }
#
#   scene = list(camera = list(eye = list(x = -1.5, y = -1.5, z = 1.5)))
#
#   tmpColors <- RColorBrewer::brewer.pal(11,"Spectral")
#   tmpColors[6] <- 'yellow'
#
#   magrittr::`%>%`(
#     magrittr::`%>%`(
#       magrittr::`%>%`(
#         plotly::plot_ly(plotData,
#                         x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
#                         type = 'scatter3d', mode = 'lines',
#                         color = ~as.factor(FDRep),
#                         colors = tmpColors),
#         plotly::layout(
#           scene = list(
#             yaxis = list(title = "EvalRange"),
#             xaxis = list(title = "FD Sims"),
#             zaxis = list(title = "Value")
#           ))),
#       plotly::layout(title = titleVal, scene = scene)),
#     plotly::layout(showlegend = FALSE))
#
#   #plotly::plot_ly(plotData,
#   #        x = ~as.factor(FDRep), y = ~evalRange, z = ~Value,
#   #        type = 'scatter3d', mode = 'lines',
#   #        color = ~as.factor(FDRep),
#   #        colors = tmpColors) %>%
#   #  #colors = c("red", "yellow", "blue")) %>%
#   #  #colors='Spectral') %>%
#   #  plotly::layout(
#   #    scene = list(
#   #      yaxis = list(title = "EvalRange"),
#   #      xaxis = list(title = "FD Sims"),
#   #      zaxis = list(title = "Value")
#   #    )) %>%
#   #  plotly::layout(title = titleVal, scene = scene) %>%
#   #  plotly::layout(showlegend = FALSE)
#   ##line = list(width = 4, color = ~as.factor(FDRep),
#   ##            colorscale = list(c(0,'#BA52ED'), c(1,'#FCB040'))))
# }
