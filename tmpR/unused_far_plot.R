
change <- 0.15
locs <- c(0.5,0.75,0.8)
lengths <- c(100,250,500)
nSims <- 1000
dat_res <- 24
alpha_val <- 0.05

results <- data.frame()

for(c in 1:length(locs)){
  tmp<- readRDS(paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\',
                 'meanFAR1_results_c',c,'.rds'))
  results <- rbind(results,tmp)
}

results <- results %>% na.omit %>%
  pivot_longer(cols = TnSim:MeanTnResid)
results100 <- results[results$length==100,]
results500 <- results[results$length==500,]

ggplot(mapping=aes(x=locs,
                   y=value,
                   group=name,
                   color=name)) +
  geom_line(data=results100,linewidth=2) +
  geom_hline(aes(yintercept=0.05),linewidth=2,
             linetype='dotted',color='gray') +
  geom_point(data=results100, size=4) +
  geom_point(data=results500, size=4) +
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        legend.justification = c(0, 0),
        legend.position = c(0.6, 0.7),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  xlab(expression(tau)) +
  ylab('Rejection Rate') #+
  # scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
  #                    labels=c('Tn','Mn','HRZ-22'),
  #                    breaks=c("TnSim","MnSim","cov"))
