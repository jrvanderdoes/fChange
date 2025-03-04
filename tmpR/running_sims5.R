
library(devtools); load_all()

change <- 0.15
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100,250,500)
nSims <- 1000
dat_res <- 24
alpha_val <- 0.05


for(c in 1:length(locs)){
  cat(paste0('\n - Change ',locs[c],' - '))
  results <- data.frame('locs'=rep(locs[c],length(lengths)),
                        'length'=lengths,
                        'TnSim'=NA,
                        'MnSim'=NA,
                        'Mean'=NA,
                        'MeanResid'=NA,
                        'MeanTn'=NA,
                        'MeanTnResid'=NA)
  for(l in c(1,3)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- readRDS(paste0(path,'Loc/',
                           'AltCase_LocChange_',locs[c],
                           '_len',lengths[l],'_data.rds'))
    tmp_results <- data.frame('TnSim'=rep(NA,nSims),
                              'MnSim'=NA,
                              'Mean'=NA,
                              'MeanResid'=NA,
                              'MeanTn'=NA,
                              'MeanTnResid'=NA)
    for(s in 1:nSims){
      cat(paste0(s,', '))

      set.seed(123*c+123^2*l+s)
      # Tn
      Tn_sim <- detect_changepoint_final_Tn(data[[s]], h=0,space='BM')

      # Mn
      Mn_sim <- detect_changepoint_final_Mn(data[[s]], h=0,space='BM')

      # Mean
      tmp_dat <- data[[s]] -
        data.frame(matrix(0,nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2),
                   matrix(changes[c],nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2))
      mean_val <- mean_change(data[[s]],inc.pval = T)
      mean_val1 <- mean_change(data[[s]],data1=tmp_dat,inc.pval = T)
      mean_valTn <- mean_change_Tn(data[[s]],inc.pval = T)
      mean_valTn1 <- mean_change_Tn(data[[s]],data1=tmp_dat,inc.pval = T)
      tmp_results[s,] <- c(Tn_sim$pval,
                           Mn_sim$pval,
                           mean_val[2],
                           mean_val1[2],
                           mean_valTn[2],
                           mean_valTn1[2])
    }

    saveRDS(tmp_results,
            paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\tmpRes_c',
                   c,'_l',l,'.rds'))

    results[l,] <- c(locs[c],lengths[l],colMeans(tmp_results<=alpha_val))
    saveRDS(results,
            paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\',
                   'meanFAR1_results_c',c,'.rds'))
  }

  saveRDS(results,paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\',
                         'loc_complete_c',c,'.rds'))
}
