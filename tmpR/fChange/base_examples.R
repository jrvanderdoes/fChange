library(devtools)
library(readxl)
library(tidyverse)
load_all()

# ### SETUP DATA (Will be saved)
# # Cancer
# WHO_bcanc <- read_excel("C:/Users/jerem/Downloads/WHO_bcanc.xlsx",
#                         sheet = "USA")
# # WHO_bcanc <- WHO_bcanc[WHO_bcanc$`Age Group`!+'[85+]',]
# WHO_bcanc_wide <- WHO_bcanc %>%
#   pivot_wider(id_cols = `Age Group`, names_from = Year,
#               values_from = `Percentage of cause-specific deaths out of total deaths`)
# cancer_funts <- funts(WHO_bcanc_wide[,-1])
#
# # Electricity
# electricity_funts <- funts(electricity,intraobs = 1:24)

# # Rates
# UScurveyields <- read.csv("C:/Users/jerem/Downloads/UScurveyields.csv")
# colnames(UScurveyields) <- c('Date','M1','M2','M3','M4','M6',
#                              'M12','M24','M36','M60','M84',
#                              'M120','M240','M360')
# UScurveyields$Date <- as.Date(UScurveyields$Date,'%m/%d/%y')
# rates <- funts(t(UScurveyields[,-1]),labels = UScurveyields[,1],
#                intraobs = c(1,2,3,4,6,12,24,36,60,84,120,240,360))
# rates_funts <- funts(rates$data[,colSums(is.na(rates$data))<nrow(rates$data)],
#                      labels = format(rates$labels[colSums(is.na(rates$data))<nrow(rates$data)],'%Y-%m'),
#                      intraobs = rates$intraobs)

cancer <- dfts(cancer)
### Summary
set.seed(1234)
png('C:/Users/jerem/Downloads/data-can/summary_electricity.png')
summary(electricity)
dev.off()

### Run generate data code
electricity$labels <- format(as.Date(electricity$labels), "%b")
rates$labels <- format(as.Date(rates$labels), "%Y")

### Plot
paper_plot <- function(X){
  val_axis_title = ''
  res_axis_title = ''
  FD_axis_title = ''
  plot_title = ''
  eye = list(x = -0.6, y = -2.5, z = 1.5)
  aspectratio = list(x = 1.75, y = 1, z = 1)

  plotData <- data.frame()
  for (i in 1:ncol(X)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$intratime,
        "FDRep" = i, #X$labels[i],
        "Value" = X$data[,i]
      )
    )
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
  plot_colors[6] <- "yellow"

  plotly::plot_ly(plotData,
                         x = ~FDRep, y = ~resolution, z = ~Value,
                         type = "scatter3d", mode = "lines",
                         color = ~ as.factor(FDRep),
                         colors = plot_colors
  ) %>%
    plotly::layout(
      scene = list(
        yaxis = list(
          title = res_axis_title,
          showticklabels = F,
          ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
          tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
        ),
        xaxis = list(
          title = FD_axis_title,
          showticklabels = T,
          ticktext=.select_n(vals=X$labels, n=5),
          tickvals=.select_n(vals=1:length(X$labels), n=5),
          tickfont = list(size = 15)
        ),
        zaxis = list(
          title = val_axis_title,
          showticklabels = F
        )
      )
    ) %>%
    plotly::layout(title = plot_title, scene = scene) %>%
    plotly::layout(showlegend = FALSE)
}
paper_plot1 <- function(X,changes=NULL){
  val_axis_title = ''
  res_axis_title = ''
  FD_axis_title = ''
  plot_title = ''
  aspectratio = c(2.5, .75, 1)


  # data <- X$data
  # curve_points <- X$intratime
  valRange <- c(
    floor(min(X$data,na.rm = T)),
    ceiling(max(X$data,na.rm = T))
  )

  name <- V1 <- value <- NULL
  data1 <- X$data
  colnames(data1) <- 1:ncol(X)
  plotData <- cbind(X$intratime,data1) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = 1+1:ncol(X)) %>%
    dplyr::mutate(name=as.numeric(name)) %>%
    dplyr::rename(resolution=V1,
                  Value=value,
                  FDRep=name)
  plotData <- plotData[order(plotData$FDRep),]

  ## Setup up Colors
  #   Rainbow for no changes, colored for changes
  plotData[["color"]] <- rep(1:ncol(X), each = nrow(X))
  if (!is.null(changes)) {
    tmp_colors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
    if (length(changes) > 9) {
      tmp_colors <- rep(tmp_colors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
    }

    changes <- unique(c(1, changes, ncol(X)))
    colors_plot <- rep(tmp_colors[1], ncol(X))
    for (i in 2:(length(changes) - 1)) {
      colors_plot[changes[i]:changes[i + 1]] <- tmp_colors[i]
    }
  } else {
    colors_plot <- RColorBrewer::brewer.pal(11, "Spectral")
    colors_plot[6] <- "yellow"
    colors_plot <- grDevices::colorRampPalette(colors_plot)(ncol(X))
  }
  plotData[["color"]] <- rep(colors_plot, each = nrow(X) )

  ## Set up Tick Labels
  #   TODO:: Test more thoroughly
  r_z <- range(plotData$Value,na.rm = T)
  r_z[2] <- r_z[2]+4
  z_range <- round(r_z, -nchar(round(max(plotData$Value,na.rm = TRUE)))+1)

  scale_info <- list(
    col = "black", arrows = FALSE, cex = 1.5,
    x = list(
      at = seq(min(X$intratime),
               max(X$intratime),
               length.out = 5
      ),
      # labels=.specify_decimal(seq(max(X$intratime),
      #                                 min(X$intratime),
      #                                 length.out=5),2)
      labels = NULL
    ),
    y = list(
      at = -seq(1, ncol(X), length.out = 7),
      labels = X$labels[seq(1, ncol(X), length.out = 7)]
    ),
    z = list(
      at = seq(z_range[1], z_range[2],
               length.out = 5
      ),
      labels = NULL
    )
  )

  color <- NULL # Check note removal
  # lattice:::wireframe
  lattice::cloud(
    x = Value ~ (resolution) * (-FDRep),
    data = plotData,
    type = "l", groups = plotData$FDRep,
    par.box = c(col = "transparent"),
    par.settings =
      list(
        axis.line = list(col = "transparent"),
        superpose.line = list(col = (colors_plot))
      ),
    # screen=list(z = 90, x = -75,y=-45),
    # trellis.par.set(list(axis.text=list(cex=2)),
    #                "axis.line",list(col=NA)),
    xlim = rev(range(X$intratime)),
    zlim = valRange,
    aspect = aspectratio,
    drape = TRUE, colorkey = FALSE,
    scales = scale_info,
    xlab = list(res_axis_title, rot = 30,cex=1.75),
    ylab = list(paste0("\n", FD_axis_title), rot = -30,cex=1.75),
    zlab = list(val_axis_title, rot = 90, just = 0.75,cex=1.75),
    main = plot_title
  )
}

plotly::save_image(paper_plot(electricity),'C:/Users/jerem/Downloads/data-can/electricity.png')
plotly::save_image(paper_plot(cancer),'C:/Users/jerem/Downloads/data-can/cancer.png')

png('C:/Users/jerem/Downloads/data-can/rates.png')
paper_plot1(rates)
dev.off()

# # 3d
# pt <- plot(dfts(cancer),val_axis_title='',plot_title = NULL)
# plotly::save_image(pt,'C:/Users/jerem/Downloads/cancer.png')
# pt <- plot(electricity,val_axis_title='',plot_title = NULL)#,
# # eye = list(x = -1.75, y = -1.75, z = 1.75))
# plotly::save_image(pt,'C:/Users/jerem/Downloads/electricity.png')
# png('C:/Users/jerem/Downloads/rates.png')
# plot(rates,type = 'fast',val_axis_title='')
# dev.off()

# 2d
png('C:/Users/jerem/Downloads/data-can/cancer_2d.png')
plot(cancer,type = 'rainbow')
dev.off()
png('C:/Users/jerem/Downloads/data-can/electricity_2d.png')
plot(electricity,type = 'rainbow')
dev.off()
# png('C:/Users/jerem/Downloads/rates_2d.png')
# .plot_stack(rates_funts)
# dev.off()

### Summarize and Prep
# impute
rates_imp <- impute(rates,method='linear')

# stationary
set.seed(1234)
stationarity_test(electricity)$pvalue
stationarity_test(cancer,critical = 'resample',blocksize = 3)$pvalue
stationarity_test(rates_imp,statistic = 'Mn')$pvalue

set.seed(1234)
kpss_test(electricity)$pvalue
kpss_test(cancer,method = 'block',h = 3)$pvalue
kpss_test(rates_imp)$pvalue

# diff
set.seed(1234)
electricity_diff <- diff(electricity)
stationarity_test(electricity_diff)$pvalue
kpss_test(electricity_diff)$pvalue

# detrend
set.seed(1234)
result <- projection_model(cancer, model = 'ets', TVE=0.99)
cancer_errs <- result$errors
kpss_test(cancer_errs)$pvalue
stationarity_test(cancer_errs)$pvalue

# White noise
set.seed(1234)
portmanteau_tests(electricity_diff,test='variety')
portmanteau_tests(cancer_errs,test='variety')
# portmanteau_tests(rates_imp,test='variety')

set.seed(12345)
png('C:/Users/jerem/Downloads/data-can/cancer_acf.png')
acf(cancer,cex.lab=1.75, cex.axis=1.1)
dev.off()
png('C:/Users/jerem/Downloads/data-can/cancer_resids_acf.png')
acf(cancer_errs,cex.lab=1.75, cex.axis=1.1)
dev.off()

###########################
set.seed(1234)
fchange(cancer,method='projmean',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]
set.seed(1234)
fchange(rates_imp,method='mean',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]
set.seed(1234)
fchange(electricity,method='robustmean',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]

set.seed(1234)
fchange(cancer,method='eigenjoint',eigen_number=3,h=2*ncol(cancer)^(1/5))[c('location','pvalue')]
set.seed(1234)
fchange(cancer,method='eigensingle',eigen_number=3,h=2*ncol(cancer)^(1/5))[c('location','pvalue')]

set.seed(1234)
fchange(rates_imp,method='trace',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]

set.seed(1234)
fchange(electricity,method='covariance',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]

set.seed(1234)
fchange(electricity,method='projdistribution',critical = 'resample',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]

set.seed(1234)
fchange(rates_imp,method='characteristic',h=2*ncol(cancer)^(1/5))[c('location','pvalue')]
###########################
paper_plot_changes <- function(X,changes){
  val_axis_title = ''
  res_axis_title = ''
  FD_axis_title = ''
  plot_title = ''
  eye = list(x = -0.6, y = -2.5, z = 1.5)
  aspectratio = list(x = 1.75, y = 1, z = 1)

  plotData <- data.frame()

  # Color Group to first CP
  for (j in 1:min(changes)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$intratime,
        "FDRep" = j,
        "Color" = 1,
        "Value" = X$data[, j]
      )
    )
  }
  # Color Group from last CP
  for (j in (max(changes) + 1):ncol(X)) {
    plotData <- rbind(
      plotData,
      data.frame(
        "resolution" = X$intratime,
        "FDRep" = j,
        "Color" = length(changes) + 1,
        "Value" = X$data[, j]
      )
    )
  }

  # Color Additional Groups
  if (length(changes) > 1) {
    for (i in 2:length(changes)) {
      for (j in (changes[i - 1] + 1):changes[i]) {
        plotData <- rbind(
          plotData,
          data.frame(
            "resolution" = X$intratime,
            "FDRep" = j,
            "Color" = i,
            "Value" = X$data[, j]
          )
        )
      }
    }
  }

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  # Get Colors
  tmpColors <- RColorBrewer::brewer.pal(min(9, max(3, length(changes) + 1)), "Set1")
  if (length(changes) > 9) {
    tmpColors <- rep(tmpColors, ceiling(c(length(changes) + 1) / 9))[1:(length(changes) + 1)]
  }

  plotly::plot_ly(plotData,
                  x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                  type = "scatter3d", mode = "lines",
                  split = ~ as.factor(FDRep),
                  color = ~ as.factor(Color),
                  colors = tmpColors
  ) %>%
    plotly::layout(
      scene = list(
        yaxis = list(
          title = res_axis_title,
          showticklabels = F,
          ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
          tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
        ),
        xaxis = list(
          title = FD_axis_title,
          showticklabels = T,
          ticktext=.select_n(vals=X$labels, n=5),
          tickvals=.select_n(vals=1:length(X$labels), n=5),
          tickfont = list(size = 15)
        ),
        zaxis = list(
          title = val_axis_title,
          showticklabels = F
        )
      )
    ) %>%
    plotly::layout(title = plot_title, scene = scene) %>%
    plotly::layout(showlegend = FALSE)
}
### Binary Segmentation - Cancer
set.seed(12345)
mean_cps <- fchange(X = cancer,method = 'mean',type='segmentation',silent.binary = FALSE)
# set.seed(12345)
# dist_cps <- change(cancer,type='segmentation',silent.binary = FALSE)

plotly::save_image(paper_plot_changes(dfts(cancer), mean_cps$location),
                   'C:/Users/jerem/Downloads/data-can/cancer_changes.png')


set.seed(12345)
mean_cps_diff <- fchange(cancer_errs,type='segmentation',silent.binary = FALSE)
set.seed(12345)
mean_cps_diff <- fchange(cancer_errs,method = 'mean',
                        type='segmentation',silent.binary = FALSE)

plotly::save_image(paper_plot(cancer_errs),
                   'C:/Users/jerem/Downloads/data-can/cancer_diff_changes.png')


### Elbow Method - Rates
set.seed(123456)
elbow_rates <- fchange(X = rates_imp,type = 'elbow',max_changes = 20)
elbow_changes <- elbow_rates$suggestion$changes
elbow_changes <- elbow_changes[order(elbow_changes)]
png('C:/Users/jerem/Downloads/data-can/rates_varplot.png')
elbow_rates$suggestion$plot
dev.off()
png('C:/Users/jerem/Downloads/data-can/rates_changes.png')
paper_plot1(rates_imp,changes=elbow_changes)
dev.off()


### Binary Segmentation - Elec
paper_interval_plot <- function(X, intervals) {
  plot_title = NULL
  val_axis_title = NULL
  res_axis_title = NULL
  FD_axis_title = NULL
  eye = list(x = -0.6, y = -2.5, z = 1.5)
  aspectratio = list(x = 1.75, y = 1, z = 1)
  highlight_changes = TRUE
  int.gradual=FALSE

  # Setup
  if(is.null(intervals)) stop('The variable `intervals` cannot be NULL.',
                              call. = FALSE)

  X <- dfts(X)
  if(!is.null(intervals)){
    changes <- unique(c(0,intervals[,1],ncol(X)))
    changes <- changes[order(changes)]

    # Get Colors
    plot_colors <- RColorBrewer::brewer.pal(
      min(9, max(length(changes) - 1,3)), "Set1")

    if (length(changes) > 9) {
      plot_colors <- rep(plot_colors,
                         ceiling(c(length(changes) + 1) / 9))[1:(length(changes) - 1)]
    }

  }else{
    changes <- 0:ncol(X)

    plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
    plot_colors[6] <- "yellow"
  }


  # Plot Parameters
  plotData <- plotData_mean <- data.frame()

  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  if(is.null(res_axis_title)){
    res_axis_title <- ''
  }
  if(is.null(FD_axis_title)){
    FD_axis_title <- ''
  }
  if(is.null(val_axis_title)){
    val_axis_title <- ''
  }

  # Color and Group (Always at least 2 changes - start/end)
  for(i in 1:(length(changes)-1)){
    region <- (changes[i]+1):changes[i+1]
    means <- rowMeans(X$data[, region, drop=FALSE])
    for(j in region){
      plotData <- rbind(
        plotData,
        data.frame(
          "resolution" = X$intratime,
          "FDRep" = j,#X$labels[j],
          "Color" = plot_colors[i],
          "Value" = X$data[, j]
        )
      )

      plotData_mean <- rbind(
        plotData_mean,
        data.frame(
          "resolution" = X$intratime,
          "FDRep" = j,#X$labels[j],
          "Color" = plot_colors[i],
          "Value" = means
        )
      )
    }

  }

  # Color Confidence Intervals and Changes
  if(!is.null(intervals)){

    if(!int.gradual){
      for(i in 1:nrow(intervals)){
        # Lower
        int <- floor(intervals[i,2]):(round(intervals[i,1]))
        plotData_mean[plotData_mean$FDRep %in% int, "Color"] <- 'gray'
        #X$labels[int], "Color"] <- 'gray'
        #colorRampPalette(c("lightgray", "darkgray"))(length(int))

        # Upper
        int <- (round(intervals[i,1])):round(intervals[i,3])
        plotData_mean[plotData_mean$FDRep %in% int, "Color"] <- 'gray'
        #X$labels[int], "Color"] <- 'gray'

        # Changes
        if(highlight_changes){
          plotData_mean[plotData_mean$FDRep %in%
                          c(round(intervals$change[i]) + 0:1),
                        # X$labels[round(intervals$change[i]) + 0:1],
                        "Color"] <- 'black'

        }
      }

      if(highlight_changes){
        plot_colors <- c(plot_colors,'black','gray')
      }else{
        plot_colors <- c(plot_colors,'gray')
      }


      return_plot <- plotly::plot_ly(plotData,
                                     x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                                     type = "scatter3d", mode = "lines",
                                     split = ~ as.factor(FDRep),
                                     color = ~ as.factor(Color),
                                     colors = plot_colors,
                                     opacity=0.15
      ) %>%
        plotly::add_lines(data=plotData_mean,
                          x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                          type = "scatter3d", mode = "lines",
                          split = ~ as.factor(FDRep),
                          color = ~ as.factor(Color),
                          colors = plot_colors,
                          opacity=1
        ) %>%
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              showticklabels = F,
              ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
              tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
            ),
            xaxis = list(
              title = FD_axis_title,
              showticklabels = T,
              ticktext=.select_n(vals=X$labels, n=5),
              tickvals=.select_n(vals=1:length(X$labels), n=5),
              tickfont = list(size = 15)
            ),
            zaxis = list(
              title = val_axis_title,
              showticklabels = F
            )
          )
        ) %>%
        plotly::layout(title = plot_title, scene = scene) %>%
        plotly::layout(showlegend = FALSE)

    }else{

      for(i in 1:nrow(intervals)){
        change <- round(intervals$change[i])
        lower <- floor(intervals$lower[i])
        upper <- floor(intervals$upper[i])

        # Lower Ramp
        # rgb <- pmin(255,grDevices::col2rgb(plot_colors[i])*1.25)/255
        # light_plot_col <- grDevices::rgb(rgb[1],rgb[2],rgb[3])
        color_ramp <-
          grDevices::colorRampPalette(#c('darkgray','darkgray'))(
            c('darkgray',plot_colors[i]))(
              change - lower+2)
        if( change >= lower ){
          for (j in change:lower) {
            # Ensure not another break/interval
            curr_color <- plotData[plotData$FDRep==j,"Color"][1]
            #X$labels[j],"Color"][1]
            use_color <- ifelse(curr_color == plot_colors[i],
                                color_ramp[change-j+1],
                                min(curr_color, color_ramp[change-j+1]))

            # Save colors
            plotData[plotData$FDRep==j, "Color"] <- use_color
            plotData_mean[plotData_mean$FDRep==j, "Color"] <- use_color
            # plotData[plotData$FDRep==X$labels[j], "Color"] <- use_color
            # plotData_mean[plotData_mean$FDRep==X$labels[j], "Color"] <- use_color
          }
        }

        # # Upper Ramp
        color_ramp <-
          grDevices::colorRampPalette(
            c('darkgray',plot_colors[i+1]))(
              upper - change+2)
        if( (change+1) <= upper ){
          for (j in (change+1):upper) {
            # Ensure not another break/interval
            curr_color <- plotData[plotData$FDRep==j,"Color"][1]
            #X$labels[j],"Color"][1]
            use_color <- ifelse(curr_color == plot_colors[i+1],
                                color_ramp[j-(change+1)+1],
                                min(curr_color, color_ramp[j-(change+1)+1]))

            # Save colors
            plotData[plotData$FDRep==j, "Color"] <- use_color
            plotData_mean[plotData_mean$FDRep==j, "Color"] <- use_color
            # plotData[plotData$FDRep==X$labels[j], "Color"] <- use_color
            # plotData_mean[plotData_mean$FDRep==X$labels[j], "Color"] <- use_color
          }
        }

        if(highlight_changes){
          plotData_mean[plotData_mean$FDRep %in% c(change+0:1),"Color"] <- 'black'
          plotData[plotData$FDRep %in% c(change+0:1),"Color"] <- 'black'
          # plotData_mean[plotData_mean$FDRep %in% X$labels[change+0:1],"Color"] <- 'black'
          # plotData[plotData$FDRep %in% X$labels[change+0:1],"Color"] <- 'black'
        }
      }


      # Plot

      pal <- plotData_mean$Color
      names(pal) <- as.factor(plotData_mean$FDRep)
      # pal <- setNames(pal, as.factor(plotData_mean$FDRep))
      pal <- pal[unique(names(pal))]

      return_plot <- plotly::plot_ly(plotData,
                                     x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                                     type = "scatter3d", mode = "lines",
                                     split = ~ as.factor(FDRep),
                                     color = ~ as.factor(FDRep),
                                     colors = pal,
                                     opacity=0.15
      ) %>%
        plotly::add_lines(data=plotData_mean,
                          x = ~ as.factor(FDRep), y = ~resolution, z = ~Value,
                          type = "scatter3d", mode = "lines",
                          split = ~ as.factor(FDRep),
                          color = ~ as.factor(FDRep),
                          colors = pal,#c(plot_colors,'black',ramp),
                          opacity=1
        ) %>%
        plotly::layout(
          scene = list(
            yaxis = list(
              title = res_axis_title,
              showticklabels = F,
              ticktext=round(seq(min(X$intratime),max(X$intratime),length.out=5),3),
              tickvals=.select_n(vals=seq(min(X$intratime),max(X$intratime),length.out=5), n=5)
            ),
            xaxis = list(
              title = FD_axis_title,
              showticklabels = T,
              ticktext=.select_n(vals=X$labels, n=5),
              tickvals=.select_n(vals=1:length(X$labels), n=5),
              tickfont = list(size = 15)
            ),
            zaxis = list(
              title = val_axis_title,
              showticklabels = F
            )
          )
        ) %>%
        plotly::layout(title = plot_title, scene = scene) %>%
        plotly::layout(showlegend = FALSE)
    }
  }

  return_plot
}

set.seed(123)
bs_result <- change(
  electricity, statistic = 'Tn', type='segmentation',silent.binary = FALSE)
interval <- confidence_interval(electricity, bs_result[[1]], alpha = 0.1)
plotly::save_image(paper_interval_plot(electricity,interval),
                   'C:/Users/jerem/Downloads/data-can/electricity_changes.png')
