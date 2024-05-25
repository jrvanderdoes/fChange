path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')


####################################
#   KL - Eigen Change Plot
####################################
library(tidyverse)

results <- data.frame()

changes <- c(1, 2, 4, 8, 16)#,3,6)
lengths <- c(100,500)
alpha_val <- 0.05
nSims <- 1000

for(i in 1:length(changes)){
  tmp <- readRDS(paste0(path,'Eigen/tmp/', 'results_c',i,'_covMethod.rds'))
  results <- rbind(results,tmp)
}
results <- results %>% na.omit %>%
  pivot_longer(cols = cov)
results$value <- 1-results$value
data_Tn <- data.frame('change'=changes,length=100,
                      TnSim=c(0.06,0.09,0.27,0.85,1)) %>%
  pivot_longer(cols=TnSim)
data_Tn1 <- data.frame('change'=changes,length=500,
                      TnSim=c(0.055,0.27,0.97,1,1)) %>%
  pivot_longer(cols=TnSim)
data_Mn <- data.frame('change'=changes,length=100,
                      MnSim=c(0.04,0.06,0.15,0.74,1)) %>%
  pivot_longer(cols=MnSim)
data_Mn1 <- data.frame('change'=changes,length=500,
                      MnSim=c(0.045,0.12,0.93,1,1)) %>%
  pivot_longer(cols=MnSim)

data_plot <- rbind(data_Tn,data_Tn1,
                   data_Mn,data_Mn1,
                   results)

results_plot <- data_plot[data_plot$change<=8,]
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]


ggplot(mapping=aes(x=change,
                   y=value,
                   group=name,
                   color=name)) +
  geom_line(data=results_plot100,linewidth=2) +
  geom_line(data=results_plot500,linewidth=2,
            linetype='dashed') +
  geom_hline(aes(yintercept=0.05),linewidth=2,
             linetype='dotted',color='gray') +
  geom_point(data=results_plot100,size=4) +
  geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.1, 0.8),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  scale_x_continuous(breaks=c(2,4,6,8),
                     limits = c(1,8)) +
  xlab(expression(Delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','HRZ-22'),
                     breaks=c("TnSim","MnSim","cov"))
