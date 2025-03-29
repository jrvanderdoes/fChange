#' Summary for dfts Object
#'
#' General summary function to view data overview. Several plots and test
#'  statistics are returned to give a general view of the data. More details
#'  can be found with more specialized functions.
#'
#' @param object A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param changes Vector of change locations, if there are any. Default is NULL.
#' @param lag.max Max lags to consider for ACF. Default is 20.
#' @param d.max Max number of dimensions for QQ-plot.
#' @param demean Boolean if data should be demeaned based on changes and create
#'  plots based on these residuals.
#' @param ... Data to pass into underlying functions like the KPSS, portmanteau,
#'  and stationary tests. In general it is recommended to not use this and
#'  instead apply the specialized functions directly.
#'
#' @return List with the elements:
#'  \enumerate{
#'    \item summary_data: summary results for the data.
#'    \item summary_plot: summary plot for the data.
#'  }
#' @export
#'
#' @seealso [base::summary()]
#'
#' @examples
#' res <- summary(electricity[,1:20], lag.max=2)
summary.dfts <- function(object, changes=NULL, lag.max=20, d.max=2, demean=FALSE, ...){
  object <- dfts(object)

  # Demean to send to correct places
  if(demean){
    if(!is.null(changes)){
      changes <- unique(c(0, changes, ncol(object)))
      object_tmp <- object

      for(d in 1:(length(changes)-1)){
        object_tmp <- dfts(object$data[,(changes[d]+1):changes[d+1]], fparam = object$fparam)
        object$data[,(changes[d]+1):changes[d+1]] <- object$data[,(changes[d]+1):changes[d+1]] - mean(object_tmp)
      }

      changes <- NULL
    }else{
      object$data <- object$data - mean(object)
    }
  }

  if(is.null(changes)) return(.summary_nochange(object, lag.max, d.max, ...))

  changes <- c(0,changes, ncol(object$data))
  changes <- changes[!is.null(changes)]
  changes <- unique(changes[order(changes)])

  #####
  ## Compute WN p-values
  p_values_wn <- data.frame(matrix(nrow=lag.max, ncol=length(changes)-1))
  data_acf <- list()

  for(cp in 1:(length(changes)-1)){
    X_cp <- dfts(object$data[,(changes[cp]+1):changes[cp+1]])
    ## WN P-values
    p_values <- rep(NA, lag.max)
    for(h1 in 1:lag.max){
      p_values[h1] <- tryCatch({
        # .single_lag_test(X_cp, lag = h)$pvalue
        .multi_lag_test(X_cp, lag = h1, ... )$pvalue
      }, error = function(e){NA} )
    }

    p_values_wn[,cp] <- p_values

    ## ACF
    data_acf[[cp]] <- acf(X_cp, lag.max=lag.max, figure=FALSE)
  }

  #####
  ## Plot white noise lags
  value <- name <- NULL

  p_values_wn_plot <-
    cbind(data.frame('lag'=1:lag.max), p_values_wn) %>%
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
  plot_qq <- .plot_distribution(object,changes = changes,d.max = d.max)

  #####
  ## Plot ACF
  SWNs <- rep(NA, length(data_acf))
  rhos <- WWNs <-
    data.frame(matrix(0, nrow=lag.max,ncol=length(data_acf)))
  for(i in 1:length(data_acf)){
    SWNs[i] <- data_acf[[i]]$SWN_bound
    max_len <- length(data_acf[[i]]$acfs)
    rhos[1:max_len,i] <- data_acf[[i]]$acfs
    WWNs[1:max_len,i] <- data_acf[[i]]$WWN_bound
  }
  rhos <- t( t(rhos) * SWNs/max(SWNs) )
  WWNs <- t( t(WWNs) * SWNs/max(SWNs) )

  . <- NULL

  rhos_long <- cbind(data.frame('lag'=1:lag.max),rhos) %>%
    tidyr::pivot_longer(cols=2:ncol(.))
  WWNs_long <- cbind(data.frame('lag'=1:lag.max),WWNs) %>%
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
  data_lines <- cbind(data.frame('Time'=object$fparam),
                      object$data) %>%
    tidyr::pivot_longer(cols = 1+1:ncol(object$data))

  if (length(changes)==2) {
    colors_plot <- as.character(1:ncol(object$data))
    for(i in 1:length(colors_plot)){
      if(nchar(colors_plot[i]) < max(nchar(colors_plot)) ){
        colors_plot[i] <- paste0(
          paste0(rep('0', max(nchar(colors_plot))-nchar(colors_plot[i])),collapse=""),
          colors_plot[i])
      }
    }
  } else{
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
    if (length(changes) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
    }

    colors_plot <- rep(tmp_colors[1], ncol(object$data))
    for (i in 2:(length(changes) - 1)) {
      colors_plot[changes[i]:changes[i + 1]] <- tmp_colors[i]
    }
  }
  data_lines$color <- rep(colors_plot,times=length(object$fparam))
  data_lines$name <- as.numeric(data_lines$name)

  Time <- color <- NULL
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
                           paste0(changes[-length(changes)],'-',changes[-1])),
               'Observations'=c(ncol(object),changes[-1]-changes[-length(changes)]),
               'kpss'=NA,
               'stationarity'=NA,
               'Resolution'=nrow(object)
               )
  for(i in 1:nrow(data_summary)){
    endpoints <- as.numeric(strsplit(data_summary$Segment[1], '-')[[1]])
    data_summary[i,'kpss'] <-
      .specify_decimal(
        kpss_test(object$data[,endpoints[1]:endpoints[2]], ...)$pvalue,
        3)
    data_summary[i,'stationarity'] <-
      .specify_decimal(
        stationarity_test(object$data[,endpoints[1]:endpoints[2]], ...)$pvalue,
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


#' Summary for dfts Object with no changes
#'
#' @inheritParams summary.dfts
#'
#' @return Plot (ggplot) summarizing the data
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .summary_nochange(generate_brownian_motion(500)
#'    \item .summary_nochange(electricity)
#'  }
#'
#' @keywords internal
#' @noRd
.summary_nochange <- function(object, lag.max = 20, d.max=2, ...){

  #####
  ## Plot white noise lags
  wn_pvalues <- rep(NA, lag.max)
  for(h in 1:lag.max){
    # wn_pvalues[h] <- .single_lag_test(object,lag = h,method = 'bootstrap')$pvalue
    wn_pvalues[h] <- .multi_lag_test(object,lag = h, ...)$pvalue
  }

  plot_whitenoise <-
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=1:lag.max, y=wn_pvalues),size=3) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=0.05),
                        linetype='dotted', col='red',linewidth=1.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=18)) +
    ggplot2::ylim(c(0,1)) +
    ggplot2::xlab('') +
    ggplot2::ylab('')

  #####
  ## Plot QQ
  plot_qq <- .plot_distribution(object,d.max=d.max)

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
  data_lines <- cbind(data.frame('Time'=object$fparam),
                      object$data) %>%
    tidyr::pivot_longer(cols = 1+1:ncol(object$data))

  colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
  colors_plot[6] <- "yellow"
  colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(object$data))
  data_lines$color <- rep(colors_plot, nrow(object$data) )
  data_lines$name <- as.numeric(data_lines$name)

  Time <- name <- value <- color <- NULL
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
                 kpss_test(object$data, ...)$pvalue,
                 3),
               'stationarity'=.specify_decimal(
                 stationarity_test(object$data, ...)$pvalue,
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
