
#' @rdname qqplot
#'
#' @export
qqplot <- function(x, ...) UseMethod("qqplot")
#' @rdname qqplot
#'
#' @export
qqplot.default <- function(x, ...) stats::qqplot(x, ...)



#' Generic Function to Compute QQ plot for funts
#'
#' Creates normal QQ plots on the principal components of functional data
#'
#' @param x funts object
#' @param TVE Numeric in [0,1] giving the total variance explained for selecting
#'  the number of principal components.
#' @param max.d Max number of principal components
#' @param qq.sep Boolean for if the QQ plots should be seperate or stacked
#' @param alpha Significant
#' @param ...
#'
#' @return ggplot2 for QQ plot
#' @export
#'
#' @examples
#' result <- qqplot(funts(electricity))
qqplot.funts <- function(x, TVE=0.95, max.d=3, qq.sep=TRUE, alpha = 0.05, ...){
  distplot(X = x, TVE = TVE, max.d = max.d, distribution = 'norm',
           qq.sep = qq.sep, alpha = alpha, ...)
}

#' Distribution plot
#'
#' @param x description
#' @param CPs description
#' @param alpha description
#' @param TVE description
#' @param max.d description
#' @param distribution description
#' @param qq.sep description
#' @param alpha description
#' @param legend description
#' @param ... description
#'
#' @return ggplot2 for plot to compare principal components to a distribution
#' @export
#'
#' @references John Fox, & Sanford Weisberg (2019). An R Companion to Applied
#'  Regression. Sage.
#'
#' @examples
#' result <- distplot(electricity)
#' result1 <- distplot(electricity, max.d=2, qq.sep=FALSE)
#' result2 <- distplot(electricity,dist='exp')
distplot <- function(X, CPs=NULL, TVE=0.95, max.d = 3,
                           distribution = "norm",
                           qq.sep = TRUE,
                           alpha = 0.05, legend = FALSE, ...){
  X <- .check_data(X)

  # Parameters
  #   All data will be of length n
  #     Since pca missing value will be fixed or already an error
  Z_alpha <- stats::qnorm(1 - alpha/2)

  # Convert X to pca
  if(is.null(CPs)){
    X_pca <- pca(X,TVE = TVE)
    D <- min(ncol(X_pca$x), max.d)

    dat <- X_pca$x[,1:D]
  }else{
    CPs <- unique(c(0, CPs, ncol(X)))
    max_n <- max(CPs[-1] - CPs[-length(CPs)])#max((CPs-lag(CPs))[-1])
    D <- length(CPs)-1
    dat <- data.frame(matrix(nrow=max_n,ncol=D))
    for(d in 1:(length(CPs)-1)){
      X_pca <- pca(funts(X$data[,(CPs[d]+1):CPs[d+1]]), TVE = 0)
      dat[,d] <- c(X_pca$x,rep(NA,length.out=max_n-length(X_pca$x)))
    }

  }


  ## Get QQ on components


  # Distribution Function
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))

  # Each components
  results <- results_coefs <- data.frame()
  for(d in 1:D){
    data_comp <- stats::na.omit(dat[,d])
    n <- length(data_comp)
    P <- stats::ppoints(n)

    ord <- order(data_comp)
    df <- data.frame('ord' = data_comp[ord],
                     'z' = q.function(P))#, ...))

    # Line
    #   c(0.25, 0.75) for x-axis on plot
    Q.x <- stats::quantile(df$ord, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75))#, ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)

    # Interval
    SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
    fit.value <- coef[1] + coef[2] * df$z
    df$upper <- fit.value + Z_alpha * SE
    df$lower <- fit.value - Z_alpha * SE

    # Label Points outside interval
    # df$label <- ifelse(df$ord > df$upper | df$ord < df$lower, X$labels[ord], "")
    results <- rbind(results, data.frame('d'=paste0('Dim-',d), df))
    results_coefs <- rbind(results_coefs,data.frame('d'=paste0('Dim-',d),'i'=coef[1],'s'=coef[2]))
  }

  if(qq.sep){
    p <- ggplot2::ggplot() +
      ggplot2::facet_grid(row=ggplot2::vars(d),scales='free_y')
  }else{
    p <- ggplot2::ggplot()
  }
  if(!legend){
    p <- p +
      ggplot2::guides(fill='none',
                      color='none')
  }

  p +
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
  # if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  # print(p)
}


