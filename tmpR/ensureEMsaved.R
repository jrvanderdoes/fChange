
elbow_method <- function(data,
                         test_function = compute_Mn,
                         trim_function = function(...){0},
                         errorType = 'L2', ... ){
  # Setup
  n <- ncol(data)
  return_data <- data.frame('CP'=NA,
                            'Var'=NA)

  # Trim & stopping criteria
  trim_amt <- trim_function(data, ...)
  nStart <- 1+trim_amt
  nEnd <- ncol(as.data.frame(data))-trim_amt
  if(nStart> nEnd) return()

  # Run First
  stats_tmp <- test_function(data, ...)
  test_stat <- stats_tmp$allValues
  return_data[1,] <- c(stats_tmp$location,
                       .compute_total_var(data, stats_tmp$location, errorType))

  # Iteratively Search
  i <- 1
  while(TRUE){
    # Setup
    i <- i + 1
    changes <- c(0,return_data$CP, n)
    changes <- changes[order(changes)]
    prev_CP <- return_data[nrow(return_data),'CP']
    prev_CP_loc <- which(changes==prev_CP)

    ## We only need to recompute for the interval changed by last CP!
    # Before
    beforePrevCP <- (changes[prev_CP_loc-1]+1):(changes[prev_CP_loc])
    test_stat[beforePrevCP] <- NA
    if(1+trim_amt < length(beforePrevCP)-trim_amt){
      fill_idx <- (1+trim_amt):(length(beforePrevCP)-trim_amt)
      stats_tmp <- test_function( data[,beforePrevCP[fill_idx]], ...)
      test_stat[beforePrevCP[fill_idx]] <- stats_tmp$allValues
    }
    # After
    afterPrevCP <- (changes[prev_CP_loc]+1):changes[prev_CP_loc+1]
    test_stat[afterPrevCP] <- NA
    if(1+trim_amt < length(afterPrevCP)-trim_amt){
      fill_idx <- (1+trim_amt):(length(afterPrevCP)-trim_amt)
      stats_tmp <- test_function( data[,afterPrevCP[fill_idx]], ...)
      test_stat[afterPrevCP[fill_idx]] <- stats_tmp$allValues
    }

    # Temp fix in case of no trim. Will update the check / var calc later
    test_stat <- ifelse(test_stat==0,NA,test_stat)

    if(sum(is.na(test_stat))==n) break

    ## Get total variance for each potential CP
    data_segments <- .split_on_NA(test_stat)

    return_data_tmp <- data.frame('CP'=NA,'Var'=NA)
    for(k in 1:length(data_segments)){
      # Find max test stat on interval
      value_max <- max(data_segments[[k]])

      # Get CP and total variance with full data
      section_max <- which(test_stat==value_max)
      tmp <- c(section_max, return_data$CP)
      tmp <- tmp[order(tmp)]

      return_data_tmp[k,] <- c(section_max,
                               .compute_total_var(data, tmp, errorType))
    }

    ## With max test-statistic on each section, take one leading to min variance
    return_data[i,] <- return_data_tmp[which.min(return_data_tmp$Var),]


    if(nrow(return_data)==(n-1)) break
  }

  # Add No Change option and compute percent change
  return_data <- rbind(c(NA, .compute_total_var(data, c(), errorType)),
                       return_data)
  return_data$Percent <- 1-return_data$Var/max(return_data$Var)

  # Define vars to remove notes
  CP <- Var <- Percent <- NULL

  var_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x=0:(length(CP)-1), y=Var)) +
    ggplot2::geom_line(ggplot2::aes(x=0:(length(CP)-1), y=Var)) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Total Variance") +
    ggplot2::theme_bw()

  per_plot <-
    ggplot2::ggplot(return_data) +
    ggplot2::geom_point(ggplot2::aes(x=0:(length(CP)-1), y=Percent)) +
    ggplot2::geom_line(ggplot2::aes(x=0:(length(CP)-1), y=Percent)) +
    ggplot2::xlab("Number of Change Points") +
    ggplot2::ylab("Percent Explained") +
    ggplot2::theme_bw()

  list('CPInfo'=return_data, 'VarPlot'=var_plot, 'PerPlot'=per_plot)
}
