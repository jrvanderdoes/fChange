path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/Loc/')

####################################
#  Loc Change Plot
####################################
alpha_val <- 0.05
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100,250,500)

for(c in 1:length(locs)){
  results <- data.frame('change'=rep(locs[c],length(lengths)),
                        'length'=lengths,
                        'TnSim'=NA,'TnBoot'=NA,'TnWelch'=NA,
                        'MnSim'=NA,'MnBoot'=NA,'Mean'=NA,'Mean_Tn'=NA)
  for(l in 1:length(lengths)){
    if(l==2) next
    tmp <- readRDS(paste0(path,'tmp/', 'loc_tmp_results_c',c,'_l',l,'.rds'))
    tmp1 <- readRDS(paste0(path,'tmp/', 'loc_tmp_results_c',c,'_l',l,'_MeanTn.rds'))
    results[l,] <- c(locs[c],lengths[l],colMeans(tmp<=alpha_val),mean(tmp1<=alpha_val))
  }
  saveRDS(results,
          paste0(path,'tmp/', 'loc_complete_c',c,'.rds'))
}

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
# data_plot$name2 <- ifelse(substr(data_plot$name,1,2)=='Me',
#                           'Mean',
#                           substr(data_plot$name,1,2))
# data_plot$name3 <- as.factor(stringr::str_replace_all(
#   stringr::str_replace_all(data_plot$name,'Tn',''),'Mn',''))

results_plot <- data_plot
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
             linetype='dashed') +
  geom_point(data=results_plot100,size=4) +
  geom_point(data=results_plot500,size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20)) +#,
  #legend.text =  element_blank(),#element_text(size=24),
  #legend.title = element_blank()) +#,
  #legend.position = c(0.87, 0.25)) +#,
  #legend.background = element_rect(fill = "white", color = "black")) +
  xlab('Location Placement') +
  ylab('Rejection Rate') #+
# guides(color="none") +
# scale_x_reverse(breaks=c(20,15,10,5,0))
