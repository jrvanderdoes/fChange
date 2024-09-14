#### PLOT RES


library(tidyverse)
load_all()

nSims <- 1000
length <- 100
dat_res <- c(20,30,40,50,60)
data <- data.frame('res'=dat_res,
                   'Tn'=NA,'TnQuant'=NA,
                   'Mn'=NA, 'MnQuant'=NA)

for(res in 1:length(dat_res)){
  results_tmp <- readRDS(
    paste0('C:\\Users\\jerem\\Downloads\\PlayHere\\res',
           res, '_100.rds'))
  data[res,-1] <- colMeans(results_tmp)
}

data_plot <- data %>%
  pivot_longer(cols = Tn:MnQuant) %>%
  mutate('resname'=paste0(res,name))
ggplot(data_plot,
       aes(x=res,y=value,color=name,fill=name,group=resname)) +
  geom_boxplot(alpha=0.25) +
  # stat_summary(geom = "errorbar",
  #              fun.min = mean, fun = mean,
  #              fun.max = mean, linewidth = 0.75,
  #              linetype='dotted') +
  theme_bw() +
  #ylim(c(0,50)) +
  scale_y_log10() +
  ggplot2::scale_color_manual(
    values=c(scales::hue_pal()(7))) +
  ggplot2::scale_fill_manual(
    values=scales::hue_pal()(7)) +
  ylab('Log Errors') +
  xlab(NULL) +
  guides(fill=guide_legend(title="Metric"),
         color=guide_legend(title="Metric"))
