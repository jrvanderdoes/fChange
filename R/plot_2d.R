#' Plot data as stack on 2d plot
#'
#' @inheritParams summary.funts
#'
#' @return Plot (ggplot2) with colored observations / regions
#'
#' @noRd
#' @keywords internal
#'
#' @examples
#' # results <- .plot_rainbow(funts(electricity))
#' # results1 <- .plot_rainbow(funts(electricity), CPs=c(50,100,200))
.plot_rainbow <- function(object, CPs=NULL){
  object <- funts(object)
  data <- object$data
  plot_data <-
    tidyr::pivot_longer(data = cbind(data.frame('Time'=object$intraobs),data),
                        cols = 1+1:ncol(data))
  # Colors
  ## Setup up Colors
  #   Rainbow for no CPs, colored for CPs
  if (!is.null(CPs)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
    if (length(CPs) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
    }
    tmp_colors <- tmp_colors[1:(length(CPs)+1)]

    CPs <- unique(c(0, CPs, ncol(data)))
    colors_plot <- rep(tmp_colors[1], ncol(data))
    for (i in 2:(length(CPs) - 1)) {
      colors_plot[(CPs[i]+1):CPs[i + 1]] <- tmp_colors[i]
    }
  } else {
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(data))
  }
  plot_data$color <- rep(colors_plot,times=nrow(data))

  # Order
  max_val <- max(nchar(plot_data$name))
  for(i in 1:nrow(plot_data)){
    if(nchar(plot_data$name[i]) < max_val ){
      plot_data$name[i] <- paste0(
        paste0(rep('0', max_val-nchar(plot_data$name[i])),collapse=""),
        plot_data$name[i])
    }
  }

  # Setup Means
  if (!is.null(CPs)) {
    means_mat <- matrix(nrow=nrow(data),ncol=length(CPs)-1)
    for(i in 1:(length(CPs)-1)){
      means_mat[,i] <- rowMeans(data[,(CPs[i]+1):CPs[i+1],drop=FALSE])
    }
    colnames(means_mat) <- unique( plot_data$color )
  } else{
    means_mat <- matrix(rowMeans(data), nrow=nrow(data), ncol=1)
  }
  means_mat <-
    tidyr::pivot_longer(data = cbind(data.frame('Time'=object$intraobs),
                                     means_mat),
                        cols = 1+1:ncol(means_mat))
  colnames(means_mat)<- c('Time','color','value')
  means_mat$name <- rep(1:(max(2,length(CPs))-1), times=nrow(data))


  ## Plot
  # For warnings
  Time <- value <- name <- color <- NULL
  ggplot2::ggplot(plot_data) +
    ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
                                    group=name,
                                    color=color),
                       alpha=0.5, linetype='dashed') +
    ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
                                    group=name,
                                    color=color),
                       data = means_mat, linewidth=1.5,
                       alpha=0.75, linetype='solid') +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::guides(color="none",alpha='none') +
    ggplot2::scale_color_manual(values = plot_data$color,
                                breaks = plot_data$color) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::xlim(range(object$intraobs))
}


#' Distribution Plot
#'
#' Plot functional data as 2-dimensional data with the mean(s) and the related
#'  distribution(s) given.
#'
#' @inheritParams summary.funts
#' @param alpha Significance
#'
#' @return Plot (ggplot2) with colored observations / regions
#'
#' @noRd
#' @keywords internal
#'
#' @examples
#' #results <- .plot_banded(funts(electricity))
#' #results1 <- .plot_banded(funts(electricity), CPs=c(50,100,200))
.plot_banded <- function(object, CPs=NULL, alpha=0.05){
  # data <- object$data
  if(!is.null(CPs)) CPs <- unique(c(0, CPs, ncol(object$data)))

  # Select subsets to show
  if (!is.null(CPs)) {
    subset_data <- data.frame(matrix(nrow=nrow(object$data),ncol=length(CPs)-1))
    for (i in 1:(length(CPs) - 1)) {
      # subsets <- c(subsets,round(stats::median((CPs_tmp[i]+1):CPs_tmp[i+1])))
      subset_data[,i] <- rowMeans(object$data[,(CPs[i]+1):CPs[i+1]])
    }
  } else{
    subset_idxs <- round(seq(1,ncol(object$data),length.out=10))
    subset_data <- object$data[,subset_idxs]
  }

  plot_data <- tidyr::pivot_longer(
    data = cbind(data.frame('Time'=object$intraobs),subset_data),
    cols = 1+1:ncol(subset_data))
  # Colors
  ## Setup up Colors
  #   Rainbow for no CPs, colored for CPs
  if (!is.null(CPs)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs)-1)), "Set1")
    if (length(CPs) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs) - 1) / 9))[1:(length(CPs)-1)]
    } else if(length(CPs)<=3){
      tmp_colors <- tmp_colors[1:(length(CPs)-1)]
    }
    colors_plot <- tmp_colors

    # CPs_tmp <- unique(c(1, CPs, ncol(data)))
    # colors_plot <- rep(tmp_colors[1], length(CPs)+1)
    # # for (i in 2:(length(CPs_tmp) - 1)) {
    # #   colors_plot[CPs_tmp[i]:CPs_tmp[i + 1]] <- tmp_colors[i]
    # for (i in 2:(length(CPs_tmp) - 1)) {
    #   colors_plot[i-1] <- tmp_colors[i]
    # }
  } else {
    # Get rainbow for all data then subset to used values
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(object$data))

    colors_plot <- colors_plot[subset_idxs]
  }
  plot_data$color <- rep(colors_plot,times=nrow(subset_data))

  # Order
  max_val <- max(nchar(plot_data$name))
  for(i in 1:nrow(plot_data)){
    if(nchar(plot_data$name[i]) < max_val ){
      plot_data$name[i] <- paste0(
        paste0(rep('0', max_val-nchar(plot_data$name[i])),collapse=""),
        plot_data$name[i])
    }
  }

  # Get Mean Bounds
  if(!is.null(CPs)){
    # warnings
    name <- group <- group1 <- group2 <- Time <- color <- type <- value <- . <- NULL

    bounds <- data.frame('Time'=object$intraobs)
    for (i in 1:(length(CPs) - 1)) {
      bounds_tmp <-
        t(apply(object$data[,(CPs[i]+1):CPs[i+1]], MARGIN = 1,
                stats::quantile, probs=c(alpha/2,1-alpha/2),na.rm=TRUE))
      colnames(bounds_tmp) <- paste0('X',i,c('_L','_U'))
      bounds <- cbind(bounds, bounds_tmp)
    }
    bounds_long <- bounds %>%
      tidyr::pivot_longer(cols = colnames(.)[-1]) %>%
      dplyr::mutate('group'=strsplit(name, '_')) %>%
      tidyr::unnest_wider(col = group,names_sep = '') %>%
      dplyr::mutate(group=group1, type=group2, group2=NULL, group1=NULL) %>%
      merge(., plot_data[,c('Time','name','color')],
            by.x = c('Time','group'), by.y = c('Time','name'),all = TRUE)
    bounds_wide <- bounds_long %>%
      tidyr::pivot_wider(id_cols = c(Time,group,color),
                         names_from = type,values_from = value)
  }else{
    bounds <- quantile.funts(object,
                             probs = c(alpha/2,1-alpha/2),
                             na.rm=TRUE)
    colnames(bounds) <- c('L','U')
    bounds_long <- cbind(data.frame('Time'=object$intraobs), bounds) %>%
      tidyr::pivot_longer(cols = colnames(.)[-1]) %>%
      dplyr::mutate('group'='X1',color='gray')
    bounds_wide <- bounds_long %>%
      tidyr::pivot_wider(id_cols = c(Time,group,color),
                         names_from = name,values_from = value)
  }

  # Plot
  # warnings
  L <- U <- NULL

  max_val <- max(bounds_long$value,plot_data$value)
  min_val <- min(bounds_long$value,plot_data$value)
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      ggplot2::aes(x=Time,
                   ymin=L, ymax=U,group=group, fill=color),
      data = bounds_wide,, alpha=0.25) +
    ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)),
                       linewidth=1.5, color='black', linetype='dashed', alpha=0.25) +
    ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
                                    group=name,
                                    color=color),
                       data = plot_data, linewidth=1) +
    ggplot2::ylim(c(min_val, max_val)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::guides(color="none",alpha='none',fill='none') +
    ggplot2::scale_color_manual(values = plot_data$color,
                                breaks = plot_data$color) +
    ggplot2::scale_fill_manual(breaks = bounds_wide$color,
                               values=bounds_wide$color) +
    ggplot2::xlab('') +
    ggplot2::ylab('')
}

