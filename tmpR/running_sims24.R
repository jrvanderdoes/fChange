## PLot EIgn
library(devtools); load_all()
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')


changes <- c(1, 2, 4, NA, 8)
lengths <- c(100,250,500)
alpha_val <- 0.05
nSims <- 1000
results100 <- results500 <-
  data.frame('changes'=changes,'Tn'=NA,'Mn'=NA, 'Cov'=NA)

for(c in c(1,2,3,5)){
  change_res <- data.frame('Tn'=rep(NA,nSims),
                           'Mn'=NA, 'Cov'=NA)
  if(c!=3){
    for(l in 1:5){
      tmp_results <-
        readRDS(paste0(path,'Eigen/100k0/',
                       'c',c,'_',l,'_100.rds'))
      change_res[1:200+(l-1)*200,] <- tmp_results[1:200+(l-1)*200,]
    }
  }else{
    change_res <-
      readRDS(paste0(path,'Eigen/100k0/',
                     'c',c,'_100.rds'))
  }
  results100[c,-1] <- c(colMeans(change_res[,1:2]<=0.05),colMeans(!is.na(change_res[3])))


  if(c !=5){
    change_res <- data.frame('Tn'=rep(NA,nSims),
                             'Mn'=NA, 'Cov'=NA)
    for(l in 1:5){
      tmp_results <-
        readRDS(paste0(path,'Eigen/500k0/',
                       'c',c,'_',l,'_500.rds'))
      change_res[1:200+(l-1)*200,] <- tmp_results[1:200+(l-1)*200,]
    }
    results500[c,-1] <- c(colMeans(change_res[,1:2]<=0.05),colMeans(!is.na(change_res[3])))
  }else{
    results500[c,-1] <- rep(1,3)
  }
}


results_plot100 <- results100 %>% na.omit %>%
  pivot_longer(cols = Tn:Cov)
results_plot500 <- results500 %>% na.omit %>%
  pivot_longer(cols = Tn:Cov)



plot_export <-
  ggplot(mapping=aes(x=changes,
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
        legend.position = c(0.7, 0.1),
        legend.text =  element_text(size=20),
        legend.title = element_blank()) +
  scale_x_continuous(breaks=c(2,4,6,8),
                     limits = c(1,8)) +
  xlab(expression(Delta)) +
  ylab('Rejection Rate') +
  scale_color_manual(values = c(scales::hue_pal()(2), 'black'),
                     labels=c('Tn','Mn','HRZ-22'),
                     breaks=c("Tn","Mn","Cov"))

png('C:\\Users\\jerem\\Downloads\\PlayHere\\cov_change.png')
print(plot_export)
dev.off()
