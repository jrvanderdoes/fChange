#' Summary for funts Object
#'
#' @param object funts object
#' @param CPs CPs, if there are any. Default is NULL.
#' @param max.lag Max lags to consider for ACF. Default is 20.
#' @param rainbow Observations shown in rainbow plot. Defaults to subset,
#'  which summarizes the values rather than showing all observations.
#'
#' @return Plot (ggplot) summarizing the data
#' @export
#'
#' @examples
#' summary(funts(electricity))
#' summary(funts(electricity), CPs=c(50,200,300))
#' summary(funts(electricity), rainbow='full')
#' summary(funts(electricity), CPs=c(50,200,300), rainbow='full')
summary.funts <- function(object, CPs=NULL, max.lag=20,
                          rainbow=c('subset','full'),
                          ...){
  object <- .check_data(object)

  if(is.null(CPs)) return(.summary_nochange(object, max.lag, rainbow=rainbow, ...))

  CPs <- c(0,CPs, ncol(object$data))
  CPs <- CPs[!is.null(CPs)]
  CPs <- unique(CPs[order(CPs)])

  #####
  ## Compute WN p-values
  p_values_wn <- data.frame(matrix(nrow=max.lag, ncol=length(CPs)-1))
  data_acf <- list()

  for(cp in 1:(length(CPs)-1)){
    X_cp <- funts(object$data[,(CPs[cp]+1):CPs[cp+1]])
    ## WN P-values
    p_values <- rep(NA, max.lag)
    for(h in 1:max.lag){
      p_values[h] <- tryCatch({
        # single_lag_test(X_cp, lag = h)$p_value
        multi_lag_test(X_cp, lag = h, method = 'iid')$p_value
      }, error = function(e){NA} )
    }

    p_values_wn[,cp] <- p_values

    ## ACF
    data_acf[[cp]] <- acf(X_cp, figure=FALSE)
  }

  #####
  ## Plot white noise lags
  p_values_wn_plot <-
    cbind(data.frame('lag'=1:max.lag), p_values_wn) %>%
    tidyr::pivot_longer(cols = colnames(p_values_wn))
  plot_whitenoise <-
    ggplot2::ggplot(p_values_wn_plot) +
      ggplot2::geom_point(ggplot2::aes(x=lag, y=value,
                                       group=name, color=name)) +
      ggplot2::geom_line(ggplot2::aes(x=lag, y=value,
                                       group=name, color=name)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),linetype='dotted', col='red') +
      ggplot2::theme_bw() +
      ggplot2::ylim(c(0,1)) +
      ggplot2::xlab('') +
      ggplot2::ylab('')  +
    ggplot2::guides(color='none')


  #####
  ## Plot rainbow
  rainbow <- match.arg(rainbow, c('subset','full'))
  if(rainbow=='full'){
    if(length(CPs)>2){
      plot_rainbow <- .plot_stack(object,CPs = CPs[-c(1,length(CPs))])
    } else{
      plot_rainbow <- .plot_stack(object)
    }
  } else if(rainbow=='subset'){
    plot_rainbow <- .plot_substack(object, CPs)
  }

  #####
  ## Plot ACF
  SWNs <- rep(NA, length(data_acf))
  rhos <- WWNs <-
    data.frame(matrix(0, nrow=max.lag,ncol=length(data_acf)))
  for(i in 1:length(data_acf)){
    SWNs[i] <- data_acf[[i]]$SWN_bound
    rhos[,i] <- data_acf[[i]]$acfs
    WWNs[,i] <- data_acf[[i]]$WWN_bound
  }
  rhos <- t( t(rhos) * SWNs/max(SWNs) )
  WWNs <- t( t(WWNs) * SWNs/max(SWNs) )

  rhos_long <- cbind(data.frame('lag'=1:max.lag),rhos) %>%
    tidyr::pivot_longer(cols=2:ncol(.))
  WWNs_long <- cbind(data.frame('lag'=1:max.lag),WWNs) %>%
    tidyr::pivot_longer(cols=2:ncol(.))

  plot_acf <-
    ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes(x=lag, y=0, yend=value, group=name, color=name),
      data=rhos_long, position=ggplot2::position_dodge(width = 0.5) ) +
    ggplot2::geom_line(
      ggplot2::aes(x=lag, y=value, group=name, color=name),
      data=WWNs_long, linetype='dashed'  ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept=SWNs,
                   color=unique(WWNs_long$name)),
      linetype='longdash') +
    ggplot2::theme_bw() +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::guides(color='none')

  #####
  ## Plot Lines
  data_lines <- cbind(data.frame('time'=object$intraobs),
                      object$data) %>%
    tidyr::pivot_longer(cols = 1+1:ncol(object$data))

  if (length(CPs)==2) {
    colors_plot <- as.character(1:ncol(object$data))
    for(i in 1:length(colors_plot)){
      if(nchar(colors_plot[i]) < max(nchar(colors_plot)) ){
        colors_plot[i] <- paste0(
          paste0(rep('0', max(nchar(colors_plot))-nchar(colors_plot[i])),collapse=""),
          colors_plot[i])
      }
    }
  } else{
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(CPs) + 1)), "Set1")
    if (length(CPs) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(CPs) + 1) / 9))[1:(length(CPs) + 1)]
    }

    colors_plot <- rep(tmp_colors[1], ncol(object$data))
    for (i in 2:(length(CPs) - 1)) {
      colors_plot[CPs[i]:CPs[i + 1]] <- tmp_colors[i]
    }
  }
  data_lines$color <- rep(colors_plot,times=length(object$intraobs))
  data_lines$name <- as.numeric(data_lines$name)

  plot_lines <- ggplot2::ggplot(data_lines,
                                ggplot2::aes(y=time,
                                    x=name,#as.Date(name)),
                                    z=value,
                                    color=as.character(color),
                                    group=name)) +
    ggplot2::theme_void() +
    stat_3D(theta=0, phi=15, geom='path') +
    ggplot2::scale_color_manual(
      breaks = unique(data_lines$color),
      values=unique(data_lines$color)) +
    ggplot2::guides(color='none')

  #####
  ## Plot Summary
  ggpubr::ggarrange(plot_lines,
                    ggpubr::ggarrange(plot_acf,
                                      plot_rainbow,
                                      ncol = 2),
                    plot_whitenoise,
                    nrow = 3 )
}


#' Summary for funts Object with no changes
#'
#' @inheritParams summary.funts
#'
#' @return Plot (ggplot) summarizing the data
#'
#' @examples
#' # .summary_nochange(
#' #   generate_brownian_motion(N = 500,
#' #      v=seq(from = 0, to = 1, length.out = 20)))
#' # .summary_nochange(funts(electricity),rainbow='full')
#' # .summary_nochange(generateFAR1(20,300))
.summary_nochange <- function(object, max.lag = 20,
                              rainbow = c('subset','full'), ...){

  #####
  ## Plot white noise lags
  wn_pvalues <- rep(NA, max.lag)
  for(h in 1:max.lag){
    # wn_pvalues[h] <- single_lag_test(object,lag = h,method = 'bootstrap')$p_value
    wn_pvalues[h] <- multi_lag_test(object,lag = h)$p_value
  }

  plot_whitenoise <-
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=1:max.lag, y=wn_pvalues)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),
                        linetype='dotted', col='red') +
    ggplot2::theme_bw() +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot rainbow TODO:: QQ-plot??
  rainbow <- match.arg(rainbow, c('subset','full'))
  if(rainbow=='full'){
    plot_rainbow <- .plot_stack(object)
  } else if(rainbow=='subset'){
    plot_rainbow <- .plot_substack(object)
  }

  #####
  ## PlotA ACF
  data_acf <- acf(object,figure=FALSE)

  plot_acf <-
    ggplot2::ggplot(
      mapping=ggplot2::aes( x=1:length(data_acf$WWN_bound) )
    ) +
    ggplot2::geom_segment( ggplot2::aes(y=0, yend=data_acf$acfs)) +
    ggplot2::geom_line(ggplot2::aes( y=data_acf$WWN_bound),
                       col='red', linetype='dashed'  ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=data_acf$SWN_bound),
                        col='blue', linetype='longdash' ) +
    ggplot2::theme_bw() +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot Lines
  data_lines <- cbind(data.frame('time'=object$intraobs),
                      object$data) %>%
    tidyr::pivot_longer(cols = 1+1:ncol(object$data))

  colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
  colors_plot[6] <- "yellow"
  colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(object$data))
  data_lines$color <- rep(colors_plot, nrow(object$data) )
  data_lines$name <- as.numeric(data_lines$name)

  plot_lines <- ggplot2::ggplot(data_lines,
                                ggplot2::aes(y=time,
                                    x=name,
                                    z=value,
                                    color=color)) +
    ggplot2::theme_void() +
    stat_3D(theta=0, phi=15, geom='path') +
    ggplot2::scale_color_manual(
      breaks = data_lines$color,
      values = data_lines$color
    ) +
    ggplot2::guides(color='none')

  #####
  ## Plot Summary
  ggpubr::ggarrange(plot_lines,
                    ggpubr::ggarrange(plot_acf,
                                      plot_rainbow,
                                      ncol = 2),
                    plot_whitenoise,
                    nrow = 3 )
}


#' Plot data as stack on 2d plot
#'
#' @inheritParams summary.funts
#'
#' @return Plot (ggplot2) with colored observations / regions
#'
#' @examples
#' # .plot_stack(funts(electricity))
#' # .plot_stack(funts(electricity), CPs=c(50,100,200))
.plot_stack <- function(object, CPs=NULL){
  object <- .check_data(object)
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
    colnames(means_mat) <- unique(plot_data$color )
  } else{
    means_mat <- matrix(rowMeans(data),nrow=nrow(data),ncol=1)
  }
  means_mat <-
    tidyr::pivot_longer(data = cbind(data.frame('Time'=object$intraobs),
                                     means_mat),
                        cols = 1+1:ncol(means_mat))
  colnames(means_mat)<- c('Time','color','value')
  means_mat$name <- rep(1:(max(2,length(CPs))-1),times=nrow(data))


  # Plot
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
    ggplot2::guides(color="none",alpha='none') +
    ggplot2::scale_color_manual(values = plot_data$color,
                                breaks = plot_data$color) +
    ggplot2::xlab('') +
    ggplot2::ylab('')
}


#' Plot data as summary stack on 2d plot
#'
#' @inheritParams summary.funts
#'
#' @return Plot (ggplot2) with colored observations / regions
#' @export
#'
#' @examples
#' # .plot_substack(funts(electricity))
#' # .plot_substack(funts(electricity), CPs=c(50,100,200))
.plot_substack <- function(object, CPs=NULL){
  # data <- object$data
  if(!is.null(CPs)) CPs <- unique(c(0, CPs, ncol(object$data)))

  # Select subsets to show
  if (!is.null(CPs)) {
    subset_data <- data.frame(matrix(nrow=nrow(object$data),ncol=length(CPs)-1))
    for (i in 1:(length(CPs) - 1)) {
      # subsets <- c(subsets,round(median((CPs_tmp[i]+1):CPs_tmp[i+1])))
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
  bounds95 <- quantile(object,probs = c(0.025,0.975))
  bounds80 <- quantile(object,probs = c(0.1,0.9))

  # if(!is.null(CPs)){
  #   # Get Mean Bounds
  #   # lowerbounds95 <- upperbounds95 <-
  #   #   lowerbounds80 <- upperbounds80 <-
  #   mean_vals <-
  #     data.frame(matrix(nrow=nrow(object$data),  ncol=length(CPs)+1))
  #   CPs_tmp <- unique(c(0, CPs, ncol(data)))
  #   for(i in 1:(length(CPs_tmp)-1)){
  #     mean_vals[,i] <-
  #       apply(object$data[,(CPs_tmp[i]+1):CPs_tmp[i+1]],
  #             MARGIN = 1, mean)
  #
  #     # lowerbounds95[,i] <-
  #     #   apply(object$data[,(CPs_tmp[i]+1):CPs_tmp[i+1]],
  #     #         MARGIN = 1, quantile,probs=0.025)
  #     # upperbounds95[,i] <-
  #     #   apply(object$data[,(CPs_tmp[i]+1):CPs_tmp[i+1]],
  #     #         MARGIN = 1, quantile,probs=0.975)
  #     #
  #     # lowerbounds80[,i] <-
  #     #   apply(object$data[,(CPs_tmp[i]+1):CPs_tmp[i+1]],
  #     #         MARGIN = 1, quantile,probs=0.1)
  #     # upperbounds80[,i] <-
  #     #   apply(object$data[,(CPs_tmp[i]+1):CPs_tmp[i+1]],
  #     #         MARGIN = 1, quantile,probs=0.9)
  #   }
  #   mean_vals <- tidyr::pivot_longer(
  #     data = cbind(data.frame('Time'=object$intraobs),mean_vals),
  #     cols = 1+1:(length(CPs)+1))
  #   # lowerbounds95 <- tidyr::pivot_longer(
  #   #   data = cbind(data.frame('Time'=object$intraobs),lowerbounds95),
  #   #   cols = 1+1:(length(CPs)+1))
  #   # upperbounds95 <- tidyr::pivot_longer(
  #   #   data = cbind(data.frame('Time'=object$intraobs),upperbounds95),
  #   #   cols = 1+1:(length(CPs)+1))
  #   # lowerbounds80 <- tidyr::pivot_longer(
  #   #   data = cbind(data.frame('Time'=object$intraobs),lowerbounds80),
  #   #   cols = 1+1:(length(CPs)+1))
  #   # upperbounds80 <- tidyr::pivot_longer(
  #   #   data = cbind(data.frame('Time'=object$intraobs),upperbounds80),
  #   #   cols = 1+1:(length(CPs)+1))
  #
  #   # Plot
  #   # return_plot <-
  #     ggplot2::ggplot() +
  #     ggplot2::geom_ribbon(ggplot2::aes(x=object$intraobs,
  #                                       ymin=bounds95[,1], ymax=bounds95[,2]),
  #                          color='gray', alpha=0.25) +
  #     ggplot2::geom_ribbon(ggplot2::aes(x=object$intraobs,
  #                                       ymin=bounds80[,1], ymax=bounds80[,2]),
  #                          color='lightgray', alpha=0.25) +
  #     ggplot2::geom_line(ggplot2::aes(x=mean_vals$Time, y=mean_vals$value,
  #                                     color=mean_vals$name, group=mean_vals$name),
  #                        linewidth=1, alpha=0.75) +
  #     ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)),
  #                        linewidth=1, color='black', alpha=0.25) +
  #     # ggplot2::geom_ribbon(ggplot2::aes(x=lowerbounds95$Time,
  #     #                                   ymin=lowerbounds95$value,
  #     #                                   ymax=upperbounds95$value,
  #     #                                   fill=lowerbounds95$name,
  #     #                                   color=lowerbounds95$name,
  #     #                                   group=lowerbounds95$name),
  #     #                      alpha=0.1) +
  #     # ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)-1.96*sd.funts(object)/sqrt(ncol(object))),
  #     #                    linewidth=1, color='black', linetype='dashed', alpha=0.25) +
  #     # ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)+1.96*sd.funts(object)/sqrt(ncol(object))),
  #     #                    linewidth=1, color='black', linetype='dashed', alpha=0.25) +
  #     # ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
  #     #                                 group=name,
  #     #                                 color=color),
  #     #                    data = plot_data, linewidth=0.5) +
  #     ggplot2::theme_bw() +
  #     ggplot2::guides(color="none",alpha='none') +
  #     ggplot2::scale_color_manual(values = plot_data$color,
  #                                 breaks = plot_data$name) +
  #     ggplot2::xlab('') +
  #     ggplot2::ylab('')
  #
  # }

  ## TODO:: Make average on for bands
  # Plot
  max_val <- max(bounds95,plot_data$value)
  min_val <- min(bounds95,plot_data$value)
  ggplot2::ggplot() +
    ggplot2::ylim(c(min_val, max_val)) +
    ggplot2::geom_ribbon(ggplot2::aes(x=object$intraobs,
                                      ymin=bounds95[,1], ymax=bounds95[,2]),
                         color='gray', alpha=0.25) +
    ggplot2::geom_ribbon(ggplot2::aes(x=object$intraobs,
                                      ymin=bounds80[,1], ymax=bounds80[,2]),
                         # ymin=min(object)$data, ymax=max(object)$data),
                         color='lightgray', alpha=0.25) +
    ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)),
                       linewidth=1, color='black', linetype='dashed', alpha=0.25) +
    # ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)-1.96*sd.funts(object)/sqrt(ncol(object))),
    #                    linewidth=1, color='black', linetype='dashed', alpha=0.25) +
    # ggplot2::geom_line(ggplot2::aes(x=object$intraobs, y=mean(object)+1.96*sd.funts(object)/sqrt(ncol(object))),
    #                    linewidth=1, color='black', linetype='dashed', alpha=0.25) +
    ggplot2::geom_line(ggplot2::aes(x=Time,y=value,
                                    group=name,
                                    color=color),
                       data = plot_data, linewidth=0.5) +
    ggplot2::theme_bw() +
    ggplot2::guides(color="none",alpha='none') +
    ggplot2::scale_color_manual(values = plot_data$color,
                                breaks = plot_data$color) +
    ggplot2::xlab('') +
    ggplot2::ylab('')
}


#' Draw 3D Geoms
#'
#' This function adds 3D geoms such as points and paths to a ggplot2 plot.
#'
#' @param theta The azimuthal rotation (degrees).
#' @param phi The colatitude rotation (degrees).
#' @param geom The geom type to use *e.g. "point", "path", or "line"*.
#' @param ... Arguments passed on to layer. Often the aesthetics
#'  like color = "red" or size = 3.
#'
#' @references Acker D (2024). gg3D: 3D perspective plots for ggplot2. R
#'                package version 0.0.0.9000.
#'
#' @export
#'
#' @examples
#' dat <- funts(electricity)
#'
#' data_lines <- cbind(data.frame('time'=dat$intraobs), dat$data) %>%
#'  tidyr::pivot_longer(cols = 1+1:ncol(dat$data))
#'
#' colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
#' colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(dat$data))
#' data_lines$color <- rep(colors_plot, nrow(dat$data) )
#' data_lines$name <- as.numeric(data_lines$name)
#'
#' ggplot2::ggplot(data_lines,
#'    ggplot2::aes(y=time, x=name, z=value, color=color)) +
#' ggplot2::theme_void() +
#' stat_3D(theta=0, phi=15, geom='path') +
#' ggplot2::scale_color_manual(
#'    breaks = data_lines$color,
#'    values = data_lines$color
#' ) +
#' ggplot2::guides(color='none')
stat_3D <- function(mapping = NULL, data = NULL, geom = "point",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
  # 3d Plot prototype
  Stat3D <- ggplot2::ggproto("Stat3D",  ggplot2::Stat,
                             setup_params = function(data, params) {

                               params$xrange <- range(data$x)
                               params$yrange <- range(data$y)
                               params$zrange <- range(data$z)

                               params
                             },
                             compute_group = function(data, scales,
                                                      theta=135, phi=60,
                                                      xrange=c(0,1),
                                                      yrange=c(0,1),
                                                      zrange=c(0,1)) {
                               data <- data %>%
                                 dplyr::mutate(
                                   x = scales::rescale(x, from=xrange, to=c(0,1)),
                                   y = scales::rescale(y, from=yrange, to=c(0,1)),
                                   z = scales::rescale(z, from=zrange, to=c(0,1)))

                               pmat <- plot3D::perspbox(z=diag(2), plot=FALSE,
                                                        theta=theta, phi=phi)

                               XY <- plot3D::trans3D(
                                 x = data$x,
                                 y = data$y,
                                 z = data$z,
                                 pmat = pmat) %>%
                                 data.frame()

                               data$x <- XY$x
                               data$y <- XY$y

                               data
                             },
                             required_aes = c("x", "y", "z")
  )

  ggplot2::layer(
    stat = Stat3D, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
