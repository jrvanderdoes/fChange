#' Plot data as stack on 2d plot
#'
#' @inheritParams summary.dfts
#'
#' @return Plot (ggplot2) with colored observations / regions
#'
#' @noRd
#' @keywords internal
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_rainbow(electricity)
#'    \item .plot_rainbow(electricity, changes=c(50,100,200))
#'  }
.plot_rainbow <- function(object, changes=NULL){
  object <- dfts(object)
  data <- object$data
  plot_data <-
    tidyr::pivot_longer(data = cbind(data.frame('Time'=object$fparam),data),
                        cols = 1+1:ncol(data))
  # Colors
  ## Setup up Colors
  #   Rainbow for no changes, colored for changes
  if (!is.null(changes)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
    if (length(changes) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
    }
    tmp_colors <- tmp_colors[1:(length(changes)+1)]

    changes <- unique(c(0, changes, ncol(data)))
    colors_plot <- rep(tmp_colors[1], ncol(data))
    for (i in 2:(length(changes) - 1)) {
      colors_plot[(changes[i]+1):changes[i + 1]] <- tmp_colors[i]
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
  if (!is.null(changes)) {
    means_mat <- matrix(nrow=nrow(data),ncol=length(changes)-1)
    for(i in 1:(length(changes)-1)){
      means_mat[,i] <- rowMeans(data[,(changes[i]+1):changes[i+1],drop=FALSE])
    }
    colnames(means_mat) <- unique( plot_data$color )
  } else{
    means_mat <- matrix(rowMeans(data), nrow=nrow(data), ncol=1)
  }
  means_mat <-
    tidyr::pivot_longer(data = cbind(data.frame('Time'=object$fparam),
                                     means_mat),
                        cols = 1+1:ncol(means_mat))
  colnames(means_mat)<- c('Time','color','value')
  means_mat$name <- rep(1:(max(2,length(changes))-1), times=nrow(data))


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
    ggplot2::xlim(range(object$fparam))
}


#' Distribution Plot
#'
#' Plot functional data as 2-dimensional data with the mean(s) and the related
#'  distribution(s) given.
#'
#' @inheritParams summary.dfts
#' @param alpha Significance
#'
#' @return Plot (ggplot2) with colored observations / regions
#'
#' @noRd
#' @keywords internal
#'
#' @examples
#' #results <- .plot_banded(electricity)
#' #results1 <- .plot_banded(electricity, changes=c(50,100,200))
.plot_banded <- function(object, changes=NULL, alpha=0.05){
  # data <- object$data
  if(!is.null(changes)) changes <- unique(c(0, changes, ncol(object$data)))

  # Select subsets to show
  if (!is.null(changes)) {
    subset_data <- data.frame(matrix(nrow=nrow(object$data),ncol=length(changes)-1))
    for (i in 1:(length(changes) - 1)) {
      # subsets <- c(subsets,round(stats::median((changes_tmp[i]+1):changes_tmp[i+1])))
      subset_data[,i] <- rowMeans(object$data[,(changes[i]+1):changes[i+1]])
    }
  } else{
    subset_idxs <- round(seq(1,ncol(object$data),length.out=10))
    subset_data <- object$data[,subset_idxs]
  }

  plot_data <- tidyr::pivot_longer(
    data = cbind(data.frame('Time'=object$fparam),subset_data),
    cols = 1+1:ncol(subset_data))
  # Colors
  ## Setup up Colors
  #   Rainbow for no changes, colored for changes
  if (!is.null(changes)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes)-1)), "Set1")
    if (length(changes) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(changes) - 1) / 9))[1:(length(changes)-1)]
    } else if(length(changes)<=3){
      tmp_colors <- tmp_colors[1:(length(changes)-1)]
    }
    colors_plot <- tmp_colors

    # changes_tmp <- unique(c(1, changes, ncol(data)))
    # colors_plot <- rep(tmp_colors[1], length(changes)+1)
    # # for (i in 2:(length(changes_tmp) - 1)) {
    # #   colors_plot[changes_tmp[i]:changes_tmp[i + 1]] <- tmp_colors[i]
    # for (i in 2:(length(changes_tmp) - 1)) {
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
  if(!is.null(changes)){
    # warnings
    name <- group <- group1 <- group2 <- Time <- color <- type <- value <- . <- NULL

    bounds <- data.frame('Time'=object$fparam)
    for (i in 1:(length(changes) - 1)) {
      bounds_tmp <-
        t(apply(object$data[,(changes[i]+1):changes[i+1]], MARGIN = 1,
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
    bounds <- quantile.dfts(object,
                             probs = c(alpha/2,1-alpha/2),
                             na.rm=TRUE)
    colnames(bounds) <- c('L','U')
    bounds_long <- cbind(data.frame('Time'=object$fparam), bounds) %>%
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
    ggplot2::geom_line(ggplot2::aes(x=object$fparam, y=mean(object)),
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

