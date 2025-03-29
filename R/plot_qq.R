#' QQ Plot Generic Function
#'
#' A generic function which by produces a qq-plot of some data. By default, it
#'  uses [stats::qqplot()].
#'
#' @param x Data to examine with a qqplot.
#' @param ... Additional parameters based on the data.
#'
#' @name qqplot
#'
#' @seealso [stats::qqplot()]
#'
#' @export
qqplot <- function(x, ...) UseMethod("qqplot")
#' @rdname qqplot
#'
#' @returns **qqplot.default**: returns results from [stats::qqplot()].
#' @export
qqplot.default <- function(x, ...) stats::qqplot(x, ...)


#' Function to Compute QQ plot for dfts Objects
#'
#' **qqplot.dfts**: Creates normal QQ plots on the principal components of functional data.
#'
#' @param x A dfts object. See [dfts()].
#' @param TVE Numeric in \[0,1\] giving the total variance explained for selecting
#'  the number of principal components.
#' @param d.max Max number of principal components. No max when NULL.
#' @param alpha Significance level, alpha in \[0,1\].
#' @param changes Vector of change points.
#' @param legend Boolean indicating if legend should be shown on plot.
#' @param ... Additional parameters based on the data.
#'
#' @returns **qqplot.dfts**: ggplot2 for QQ plot.
#' @export
#'
#' @rdname qqplot
#'
#' @examples
#' result <- qqplot(electricity, d.max=3)
qqplot.dfts <- function(x, TVE=0.95, d.max=NULL, alpha = 0.05,
                         changes=NULL, legend = FALSE, ...){
  .plot_distribution(X = x, TVE = TVE, d.max = d.max, distribution = 'norm',
                     alpha = alpha, changes=changes, legend=legend, ...)
}


#' Distribution plot
#'
#' @param X A dfts object or data which can be automatically converted to that
#'  format. See [dfts()].
#' @param changes Vector of change points
#' @param TVE Total variance explained. Value in \[0,1\]
#' @param d.max Numberic for max number of pca components to examine
#' @param distribution String for distribution to consider. Any distribution with
#'  d\code{distribution} and q\code{distribution} defined can be used
#' @param alpha Significance level, alpha in \[0,1\]
#' @param legend Boolean indicating if legend should be shown
#' @param ... Additional parameters for the distribution
#'
#' @return ggplot2 for plot to compare principal components to a distribution
#'
#' @references John Fox, & Sanford Weisberg (2019). An R Companion to Applied
#'  Regression. Sage.
#'
#' @details The following examples may be useful if this (internal) function
#'  is investigated.
#'  \itemize{
#'    \item .plot_distribution(electricity)
#'    \item .plot_distribution(electricity,distribution='exp', changes=c(50,175))
#'  }
#'
#' @noRd
#' @keywords internal
.plot_distribution <- function(X, changes=NULL, TVE=0.95, d.max = 3,
                               distribution = "norm",
                               alpha = 0.05, legend = FALSE, ...){
  X <- dfts(X)

  # Parameters
  #   All data will be of length n
  #     Since pca missing value will be fixed or already an error
  Z_alpha <- stats::qnorm(1 - alpha/2)

  # Convert X to pca
  if(is.null(changes)){
    X_pca <- pca(X,TVE = TVE)
    D <- min(ncol(X_pca$x), d.max)

    dat <- X_pca$x[,1:D]
  }else{
    changes <- unique(c(0, changes, ncol(X)))
    max_n <- max(changes[-1] - changes[-length(changes)])#max((changes-lag(changes))[-1])
    D <- min(nrow(X), d.max)
    dat <- list()
    for(d in 1:(length(changes)-1)){
      X_pca <- pca(dfts(X$data[,(changes[d]+1):changes[d+1]]), TVE = TVE)
      dat[[d]] <- X_pca$x[,1:D]
      # dat[,d] <- c(X_pca$x,rep(NA,length.out=max_n-length(X_pca$x)))
    }

  }


  ## Get QQ on components


  # Distribution Function
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))

  # Each components
  results <- results_coefs <- data.frame()
  for(d in 1:D){
    n <- c()
    P <- list()

    if(!is.null(changes)){
      df <- data.frame()
      for(i in 1:(length(changes)-1)){
        data_comp <- stats::na.omit(dat[[i]][,d])
        n <- c(n,length(data_comp))
        P[[i]] <- stats::ppoints(n[i])

        ord <- order(data_comp)
        df <- rbind(df,
                    data.frame(
                      'CP' = paste0(changes[i]+1,'-',changes[i+1]),
                      'ord' = data_comp[ord],
                      'z' = q.function(P[[i]], ...)))
      }
    }else{
      data_comp <- stats::na.omit(dat[,d])
      n <- length(data_comp)
      P[[1]] <- stats::ppoints(n)

      ord <- order(data_comp)
      df <- data.frame('CP'=paste0('1-',ncol(X)),
                       'ord' = data_comp[ord],
                       'z' = q.function(P[[1]], ...))
    }

    # Line
    #   c(0.25, 0.75) for x-axis on plot
    groups <- unique(df$CP)
    for(i in 1:length(groups)){
      Q.x <- stats::quantile(df[df$CP==groups[i],'ord'], c(0.25, 0.75))
      Q.z <- q.function(c(0.25, 0.75), ...)
      b <- diff(Q.x)/diff(Q.z)
      coef <- c(Q.x[1] - b * Q.z[1], b)

      # Interval
      SE <- (coef[2]/d.function(df[df$CP==groups[i],'z'], ...)) * sqrt(P[[i]] * (1 - P[[i]])/n[i])
      fit.value <- coef[1] + coef[2] * df[df$CP==groups[i],'z']
      df[df$CP==groups[i],'upper'] <- fit.value + Z_alpha * SE
      df[df$CP==groups[i],'lower'] <- fit.value - Z_alpha * SE

      # Setup For Plot
      results_coefs <- rbind(results_coefs,
                             data.frame('CP'=groups[i],'d'=paste0('Dim-',d),
                                        'i'=coef[1],'s'=coef[2]))
    }

    # Label Points outside interval
    # df$label <- ifelse(df$ord > df$upper | df$ord < df$lower, X$labels[ord], "")
    results <- rbind(results, data.frame('d'=paste0('Dim-',d), df))
  }

  if(!is.null(changes)){
    p <- ggplot2::ggplot() +
      ggplot2::facet_grid(row=ggplot2::vars(d),scales='free_y')
  }else{
    p <- ggplot2::ggplot()
  }
  if(!legend){
    p <- p +
      ggplot2::guides(fill='none',
                      color='none')
  }else{
    if(!is.null(changes)){
      # Order legend
      p <- p +
        ggplot2::scale_fill_discrete(breaks=unique(results$CP),
                                     guide=ggplot2::guide_legend(override.aes = list(linetype = 0))) +
        ggplot2::scale_color_discrete(breaks=unique(results$CP))
    }else{
      # Order legend
      p <- p +
        ggplot2::scale_fill_discrete(breaks=unique(results$d),
                                     guide=ggplot2::guide_legend(override.aes = list(linetype = 0))) +
        ggplot2::scale_color_discrete(breaks=unique(results$d))
    }
  }

  # Remove warnings
  z <- lower <- upper <- i <- s <- CP <- NULL

  if(!is.null(changes)){
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(x=z, ymin = lower, ymax = upper, color=CP, fill=CP),
                           data=results, alpha=0.2) +
      ggplot2::geom_abline(ggplot2::aes(intercept = i, slope = s,
                                        color=CP), data= results_coefs) +
      ggplot2::geom_point(ggplot2::aes(x=z,y=ord,color=CP,fill=CP),data=results) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                     axis.text.y = ggplot2::element_blank()) +
      ggplot2::labs(fill='',color='') +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::scale_x_continuous(breaks=floor(min(results$z)):ceiling(max(results$z)))
  } else{
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(x=z, ymin = lower, ymax = upper, color=d, fill=d),
                           data=results, alpha=0.2) +
      ggplot2::geom_abline(ggplot2::aes(intercept = i, slope = s,
                                        color=d), data= results_coefs) +
      ggplot2::geom_point(ggplot2::aes(x=z,y=ord,color=d,fill=d),data=results) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = ggplot2::element_text(size=18),
                     axis.text.y = ggplot2::element_blank()) +
      ggplot2::labs(fill='',color='') +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::scale_x_continuous(breaks=floor(min(results$z)):ceiling(max(results$z)))
  }

  p
}


