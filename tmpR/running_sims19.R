path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')

####################################
#   KL - Mean Change
####################################
library(tidyverse)

results <- data.frame()
changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
lengths <- c(100,250,500)

for(i in 1:length(changes)){
  tmp <- readRDS(paste0(path,'Mean/', 'mean_complete_c',i,'.rds'))
  tmp1 <- readRDS(paste0(path,'Mean/', 'mean_complete_c',i,'_MeanMethod.rds'))

  results <- rbind(results,merge(tmp,tmp1))
}

data_plot <- results %>% na.omit %>%
  pivot_longer(cols = TnSim:Mean)
data_plot$name <- ifelse(data_plot$name=='Mean',
                         'Mean',
                         paste0(substr(data_plot$name,1,2),' ',
                                substr(data_plot$name,3,10)))

results_plot <- data_plot
results_plot <- results_plot[!(results_plot$name%in%c('Tn Boot','Mn Boot','Tn Welch')),]
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]

plot_export <- ggplot(mapping=aes(x=change,
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
        legend.position = c(0.68, 0.16),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','ARS-18'),
                     breaks=c("Tn Sim","Mn Sim","Mean"))

png('C:\\Users\\jerem\\Downloads\\mean_change.png')
print(plot_export)
dev.off()

####################################
#   KL - Eigen Change
####################################
library(tidyverse)

results <- data.frame()
changes <- c(1, 2, 4, 8, 16, 3, 6)
lengths <- c(100,250,500)

for(i in 1:length(changes)){
  tmp <- readRDS(paste0(path,'Eigen/', 'mean_complete_c',i,'.rds'))
  #tmp1 <- readRDS(paste0(path,'Eigen/', 'mean_complete_c',i,'_MeanMethod.rds'))

  results <- rbind(results,merge(tmp,tmp1))
}

data_plot <- results %>% na.omit %>%
  pivot_longer(cols = TnSim:Mean)
data_plot$name <- ifelse(data_plot$name=='Mean',
                         'Mean',
                         paste0(substr(data_plot$name,1,2),' ',
                                substr(data_plot$name,3,10)))

results_plot <- data_plot
results_plot <- results_plot[!(results_plot$name%in%c('Tn Boot','Mn Boot','Tn Welch')),]
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]


plot_export <- ggplot(mapping=aes(x=change,
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
        legend.position = c(0.7, 0.16),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(Delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','Mean'),
                     breaks=c("Tn Sim","Mn Sim","Mean"))

png('C:\\Users\\jerem\\Downloads\\cov_change.png')
print(plot_export)
dev.off()

####################################
#   KL - Dist Change
####################################
library(tidyverse)

changes <- c(20,10,5,4,3,2,1)
lengths <- c(100,250,500)
results <- data.frame()

null_path <- paste0(getwd(),'/paperFiles/Null/KL/Results/')
# Null Case
results <-
  rbind(results,
        data.frame('change'=Inf,
                   'length'=100,
                   'TnSim' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_len1.rds'))<=0.05)[2],
                   'TnBoot' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_boot_len1.rds'))<=0.05)[2],
                   'TnWelch' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_welch_len1.rds'))>=0.95)[2],
                   'MnSim' = colMeans(readRDS(
                     paste0(null_path,
                            'Mn_NullCase_len1.rds'))<=0.05)[2],
                   'MnBoot' = colMeans(readRDS(
                     paste0(null_path,
                            'Mn_NullCase_boot_len1.rds'))<=0.05)[2],
                   'Mean' = mean(readRDS(
                     paste0(null_path,
                            'meanCase_len1.rds'))<=0.05))
  )
results <-
  rbind(results,
        data.frame('change'=Inf,
                   'length'=500,
                   'TnSim' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_len3.rds'))<=0.05)[2],
                   'TnBoot' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_boot_len3.rds'))<=0.05)[2],
                   'TnWelch' = colMeans(readRDS(
                     paste0(null_path,
                            'Tn_NullCase_welch_len3.rds'))>=0.95)[2],
                   'MnSim' = colMeans(readRDS(
                     paste0(null_path,
                            'Mn_NullCase_len3.rds'))<=0.05)[2],
                   'MnBoot' = colMeans(readRDS(
                     paste0(null_path,
                            'Mn_NullCase_boot_len3.rds'))<=0.05)[2],
                   'Mean' = mean(readRDS(
                     paste0(null_path,
                            'meanCase_len3.rds'))<=0.05))

  )

for(i in 2:length(changes)){
  tmp <- readRDS(paste0(path,'Dist/', 'dist_complete_c',i,'.rds'))
  #tmp1 <- readRDS(paste0(path,'Dist/', 'mean_complete_c',i,'_MeanMethod.rds'))
  tmp$change <- changes[i]

  results <- rbind(results,data.frame(tmp,'Mean'=NA))#merge(tmp,tmp1))
}

data_plot <- results %>%
  mutate(Mean=NULL) %>%
  na.omit %>%
  pivot_longer(cols = TnSim:MnBoot)#Mean)
data_plot$name <- ifelse(data_plot$name=='Mean',
                         'Mean',
                         paste0(substr(data_plot$name,1,2),' ',
                                substr(data_plot$name,3,10)))

results_plot <- data_plot
results_plot <- results_plot[!(results_plot$name%in%c('Tn Boot','Mn Boot','Tn Welch')),]
results_plot100 <- results_plot[results_plot$length==100,]
results_plot500 <- results_plot[results_plot$length==500,]


plot_export <- ggplot(mapping=aes(x=change,
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
  xlab(expression(nu)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','Mean'),
                     breaks=c("Tn Sim","Mn Sim","Mean")) +
  scale_x_reverse(breaks=c(15,10, 8, 6, 4, 2, 1),
                  labels=c('...',10,8,6,4,2,1),
                  lim =c(15,1))

png('C:\\Users\\jerem\\Downloads\\dist_change.png')
print(plot_export)
dev.off()

#####################################################################
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/Loc/')

####################################
#  Loc Change Plot
####################################
library(tidyverse)

results <- data.frame()
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100,250,500)

for(i in 1:length(locs)){
  tmp <- readRDS(paste0(path,'loc_complete_c',i,'.rds'))
  results <- rbind(results,tmp)
}

data_plot <- results %>% na.omit %>%
  pivot_longer(cols = TnSim:Mean_Tn)
data_plot$name <- ifelse(data_plot$name=='Mean',
                         'Mean',
                         ifelse(data_plot$name=='Mean_Tn',
                                'Mean Tn',
                                paste0(substr(data_plot$name,1,2),' ',
                                       substr(data_plot$name,3,10))))

results_plot <- data_plot
results_plot <- results_plot[!(results_plot$name%in%c('Tn Boot','Mn Boot','Tn Welch')),]
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
             linetype='dotted',color='gray') +
  geom_point(data=results_plot100,size=4) +
  geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.67, 0.75),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(tau)) +
  ylab('Rejection Rate') #+
  # scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
  #                    labels=c('Tn','Mn','Mean'),
  #                    breaks=c("Tn Sim","Mn Sim","Mean"))

png('C:\\Users\\jerem\\Downloads\\loc_change.png')
print(plot_export)
dev.off()
