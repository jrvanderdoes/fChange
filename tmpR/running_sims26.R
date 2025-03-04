
####################################
#   OU - Mean Change Plot
####################################
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/OU/')


library(tidyverse)

results <- data.frame()
changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,0.4,0.5)
lengths <- c(100,250,500)

for(i in c(1,3,6,8)){
  tmp <- readRDS(paste0(path,'Mean1/', 'complete_results_c',i,'.rds'))
  colnames(tmp) <- c('change','length','Tn','Mn','Mean')
  results <- rbind(results,tmp)
}

data_plot <- results %>% na.omit %>%
  pivot_longer(cols = Tn:Mean)

results_plot <- data_plot
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]


plot_export <-
  ggplot(mapping=aes(x=change,
                     y=value,
                     group=name,
                     color=name)) +
  geom_line(data=results_plot100,linewidth=2) +
  geom_line(data=results_plot500,linewidth=2,
            linetype='dashed') +
  geom_hline(aes(yintercept=0.05),linewidth=2,
             color='gray',linetype='dotted') +
  geom_point(data=results_plot100,size=4) +
  geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.1, 0.7),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black',scales::hue_pal()(4)[4]),
                     labels=c('Tn','Mn',"ARS-18"),
                     breaks=c("Tn","Mn","Mean"))
plot_export

png('C:\\Users\\jerem\\Downloads\\PlayHere\\ou_change.png')
print(plot_export)
dev.off()



####################################
#   FAR(1) - Mean Change Plot
####################################
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/FAR1/')


library(tidyverse)

results <- data.frame()
changes <- c(0, 0.3, 0.4, 0.35, 0.25, 0.15)
lengths <- c(100,250,500)

for(i in 1:5){
  tmp <- readRDS(paste0(path,'Mean1/', 'far1_complete_c',i,'.rds'))
  results <- rbind(results,tmp)
}

data_plot <- results %>% na.omit %>%
  pivot_longer(cols = Tn:Mean)

results_plot <- data_plot
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]


plot_export <-
  ggplot(mapping=aes(x=changes,
                   y=value,
                   group=name,
                   color=name)) +
  geom_line(data=results_plot100,linewidth=2) +
  geom_line(data=results_plot500,linewidth=2,
            linetype='dashed') +
  geom_hline(aes(yintercept=0.05),linewidth=2,
             color='gray',linetype='dotted') +
  geom_point(data=results_plot100,size=4) +
  geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.1, 0.7),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black',scales::hue_pal()(4)[4]),
                     labels=c('Tn','Mn',"ARS-18"),
                     breaks=c("Tn","Mn","Mean"))
png('C:\\Users\\jerem\\Downloads\\PlayHere\\far1_change.png')
print(plot_export)
dev.off()

