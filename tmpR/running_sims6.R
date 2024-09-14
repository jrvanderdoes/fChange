library(tidyverse)
####################################
#   OU(1) - Mean Change
####################################
library(devtools); load_all()

changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
lengths <- c(100,500)
nSims <- 1000
dat_res <- 24
results <-
  data.frame('changes'=changes,'TnSim'=NA,'TnBoot'=NA,'TnWelch'=NA,
             'MnSim'=NA,'MnBoot'=NA, 'Mean'=NA)

for(c in 1:length(changes)){
  change_res <- data.frame('TnSim'=rep(NA,nSims),'TnBoot'=NA,'TnWelch'=NA,
                           'MnSim'=NA,'MnBoot'=NA, 'Mean'=NA)
  for(i in 1:5){
    tmp <- readRDS(paste0('C:/Users/jerem/Downloads/PlayHere/OU/',
                   'meanOU_tmp_results_c',c,'_l1_',i,'.rds'))
    change_res[1:200+(i-1)*200,] <- tmp[1:200+(i-1)*200,]
  }
  results[c,-1] <- colMeans(change_res<=0.05)
}

results_plot <- results %>% na.omit %>%
  mutate(TnBoot=NULL,TnWelch=NULL,MnBoot=NULL) %>%
  pivot_longer(cols = TnSim:Mean)



plot_export <-
  ggplot(mapping=aes(x=changes,
                     y=value,
                     group=name,
                     color=name)) +
  geom_line(data=results_plot,linewidth=2) +
  # geom_line(data=results_plot500,linewidth=2,
  #           linetype='dashed') +
  geom_hline(aes(yintercept=0.05),linewidth=2,
             linetype='dotted',color='gray') +
  geom_point(data=results_plot,size=4) +
  # geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.1, 0.7),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  # scale_x_continuous(breaks=c(2,4,6,8),
  #                    limits = c(1,8)) +
  xlab(expression(delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','ARS-18'),
                     breaks=c("TnSim","MnSim","Mean"))

png('C:\\Users\\jerem\\Downloads\\PlayHere\\OU_change.png')
print(plot_export)
dev.off()
