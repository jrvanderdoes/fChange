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
#' summary(funts(electricity[,1:20]), max.lag=3)
#' # summary(funts(electricity[,1:50]), CPs=c(20,38), max.lag=3)
#' # summary(funts(electricity), CPs=c(50,200,300), max.lag=5, demean = TRUE)
summary.funts <- function(object, CPs=NULL, max.lag=20, max.d=2, demean=FALSE, ...){
  object <- .check_data(object)

  # Demean to send to correct places
  if(demean){
    if(!is.null(CPs)){
      CPs <- unique(c(0, CPs, ncol(object)))
      object_tmp <- object

      for(d in 1:(length(CPs)-1)){
        object_tmp <- funts(object$data[,(CPs[d]+1):CPs[d+1]], intraobs = object$intraobs)
        object$data[,(CPs[d]+1):CPs[d+1]] <- object$data[,(CPs[d]+1):CPs[d+1]] - mean(object_tmp)
      }

      CPs <- NULL
    }else{
      object$data <- object$data - mean(object)
    }
  }

  if(is.null(CPs)) return(.summary_nochange(object, max.lag, max.d, ...))

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
    data_acf[[cp]] <- acf(X_cp, lag.max=max.lag, figure=FALSE)
  }

  #####
  ## Plot white noise lags
  p_values_wn_plot <-
    cbind(data.frame('lag'=1:max.lag), p_values_wn) %>%
    tidyr::pivot_longer(cols = colnames(p_values_wn)) %>%
    stats::na.omit()
  plot_whitenoise <-
    ggplot2::ggplot(p_values_wn_plot) +
      ggplot2::geom_point(ggplot2::aes(x=lag, y=value,
                                       group=name, color=name)) +
      ggplot2::geom_line(ggplot2::aes(x=lag, y=value,
                                       group=name, color=name)) +
      ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),linetype='dotted', col='red') +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
      ggplot2::ylim(c(0,1)) +
      ggplot2::xlab('') +
      ggplot2::ylab('')  +
    ggplot2::guides(color='none')


  #####
  ## Plot QQ
  plot_qq <- distplot(object,CPs = CPs,max.d = max.d,qq.sep = FALSE)

  #####
  ## Plot ACF
  SWNs <- rep(NA, length(data_acf))
  rhos <- WWNs <-
    data.frame(matrix(0, nrow=max.lag,ncol=length(data_acf)))
  for(i in 1:length(data_acf)){
    SWNs[i] <- data_acf[[i]]$SWN_bound
    max_len <- length(data_acf[[i]]$acfs)
    rhos[1:max_len,i] <- data_acf[[i]]$acfs
    WWNs[1:max_len,i] <- data_acf[[i]]$WWN_bound
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
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::guides(color='none')

  #####
  ## Plot Lines
  data_lines <- cbind(data.frame('Time'=object$intraobs),
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
                                ggplot2::aes(y=Time,
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
  ## Print Descriptive Statistics
  data_summary <-
    data.frame('Segment'=c(paste0('1-',ncol(object)),
                           paste0(CPs[-length(CPs)],'-',CPs[-1])),
               'Observations'=c(ncol(object),CPs[-1]-CPs[-length(CPs)]),
               'kpss'=NA,
               'stationarity'=NA,
               'Resolution'=nrow(object)
               )
  for(i in 1:nrow(data_summary)){
    endpoints <- as.numeric(stringr::str_split(data_summary$Segment[1],'-')[[1]])
    data_summary[i,'kpss'] <-
      .specify_decimal(
        compute_kpss(object$data[,endpoints[1]:endpoints[2]],
                     method="BS", boot_method='overlapping')$pvalue,
        3)
    data_summary[i,'stationarity'] <-
      .specify_decimal(
        compute_stationarity_test(object$data[,endpoints[1]:endpoints[2]])$pvalue,
        3)
  }

  #####
  ## Plot Summary
  plot_summary <- ggpubr::ggarrange(plot_lines,
                                    ggpubr::ggarrange(plot_acf,
                                                      plot_qq,
                                                      ncol = 2),
                                    plot_whitenoise,
                                    nrow = 3 )

  list('summary_data'=data_summary,
       'summary_figure'=plot_summary)
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
#' # .summary_nochange(generate_far1(20,300))
#'
#' @keywords internal
#' @noRd
.summary_nochange <- function(object, max.lag = 20, max.d=2, ...){

  #####
  ## Plot white noise lags
  wn_pvalues <- rep(NA, max.lag)
  for(h in 1:max.lag){
    # wn_pvalues[h] <- single_lag_test(object,lag = h,method = 'bootstrap')$p_value
    wn_pvalues[h] <- multi_lag_test(object,lag = h)$p_value
  }

  plot_whitenoise <-
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=1:max.lag, y=wn_pvalues),size=3) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),
                        linetype='dotted', col='red',linewidth=1.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot QQ
  plot_qq <- distplot(object,max.d=max.d,qq.sep = FALSE)

  #####
  ## Plot ACF
  data_acf <- acf(object,figure=FALSE)

  plot_acf <-
    ggplot2::ggplot(
      mapping=ggplot2::aes( x=1:length(data_acf$WWN_bound) )
    ) +
    ggplot2::geom_segment( ggplot2::aes(y=0, yend=data_acf$acfs),
                           linewidth=2) +
    ggplot2::geom_line(ggplot2::aes( y=data_acf$WWN_bound),
                       col='red', linetype='dashed',
                       linewidth=2  ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=data_acf$SWN_bound),
                        col='blue', linetype='longdash',
                        linewidth=2) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot Lines
  data_lines <- cbind(data.frame('Time'=object$intraobs),
                      object$data) %>%
    tidyr::pivot_longer(cols = 1+1:ncol(object$data))

  colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
  colors_plot[6] <- "yellow"
  colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(object$data))
  data_lines$color <- rep(colors_plot, nrow(object$data) )
  data_lines$name <- as.numeric(data_lines$name)

  plot_lines <- ggplot2::ggplot(data_lines,
                                ggplot2::aes(y=Time,
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
  ## Print Descriptive Statistics
  data_summary <-
    data.frame('Segment'=paste0('1-',ncol(object)),
               'Observations'=ncol(object),
               'kpss'=.specify_decimal(
                 compute_kpss(object$data, method="BS", boot_method='overlapping')$pvalue,
                 3),
               'stationarity'=.specify_decimal(
                 compute_stationarity_test(object$data)$pvalue,
                 3),
               'Resolution'=nrow(object)
    )

  #####
  ## Plot Summary
  plot_summary <- ggpubr::ggarrange(plot_lines,
                                    ggpubr::ggarrange(plot_acf,
                                                      plot_qq,
                                                      ncol = 2),
                                    plot_whitenoise,
                                    nrow = 3 ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0,0.2,0,0.2, "cm"))

  list('summary_data'=data_summary,
       'summary_figure'=plot_summary)
}


#' Plot data as stack on 2d plot
#'
#' @inheritParams summary.funts
#'
#' @return Plot (ggplot2) with colored observations / regions
#' @export
#'
#' @examples
#' rainbow_plot(funts(electricity))
#' rainbow_plot(funts(electricity), CPs=c(50,100,200))
rainbow_plot <- function(object, CPs=NULL){
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
#'
#' @return Plot (ggplot2) with colored observations / regions
#' @export
#'
#' @examples
#' distribution_plot(funts(electricity))
#' distribution_plot(funts(electricity), CPs=c(50,100,200))
distribution_plot <- function(object, CPs=NULL, alpha=0.05){
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
      dplyr::mutate('group'=stringr::str_split(name,'_')) %>%
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
#' data_lines <- cbind(data.frame('Time'=dat$intraobs), dat$data) %>%
#'  tidyr::pivot_longer(cols = 1+1:ncol(dat$data))
#'
#' colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
#' colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(dat$data))
#' data_lines$color <- rep(colors_plot, nrow(dat$data) )
#' data_lines$name <- as.numeric(data_lines$name)
#'
#' result <- ggplot2::ggplot(data_lines,
#'    ggplot2::aes(y=Time, x=name, z=value, color=color)) +
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
