## Temperature
library(fda)
library(tidyverse)

temperature <- fChange::Australian_Temp$Sydney
temperature <-
  temperature[,c('Year','Month','Day',
                 'Minimum.temperature..Degree.C.')]
colnames(temperature) <- c('Year','Month','Day','Temp')
temperature <- temperature %>%
  #mutate('Year'=paste0('Y',Year)) %>%
  pivot_wider(id_cols = c(Month,Day),
              names_from = Year,
              values_from = Temp) %>%
  arrange(Month,Day) %>%
  #mutate('DayOfYear'=1:nrow(.)) %>%
  #mutate('PercentofYear'=DayOfYear/max(DayOfYear)) %>%
  #relocate(PercentofYear) %>%
  #dplyr::select(-DayOfYear) %>%
  filter(!(Month==2 & Day == 29)) %>%
  dplyr::select(-c(Month,Day)) %>%
  dplyr::select(-c(`2012`)) %>%
  as.matrix() %>%
  linear_imputatation(.)

temp_bs <- complete_binary_segmentation(
  data=temperature,
  # test_statistic_function = compute_Tn,
  # changepoint_function = detect_changepoint_singleCov,
  trim_function=function(data, ...){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)},
  changepoint_function=function(data, alpha, ...){
    func_val <- generalized_resampling(X=data,
                                       blockSize=1,
                                       fn=compute_Tn, iters=1000, replace=F, ...)

    ifelse(func_val$pval < alpha,
           compute_Mn(data, which.Mn=TRUE, ...),
           NA)
  },
  M=1000,
  final_verify = T,
  silent = F,
  space='BM', h=0)

# Paper Plotting
plot_fd(temperature, showticklabels = F,
        FD_axis_title = '',val_axis_title = '',
        eval_axis_title = '',
        eye = list(x = -.25, y = -1.5, z = 0.75),
        aspectratio=list(x=1.5,y=0.5,z=0.75))
plot_fd(temperature,temp_bs, showticklabels = F,
        FD_axis_title = '',val_axis_title = '',
        eval_axis_title = '',
        eye = list(x = -.25, y = -1.5, z = 0.5),
        aspectratio=list(x=1.5,y=0.5,z=0.75))
