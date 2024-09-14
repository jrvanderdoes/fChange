path <- getwd()
path <- paste0('C:\\Users\\jerem\\Downloads\\PlayHere\\')

####################################
#   Sin Mean Change
####################################
library(devtools); load_all()

Ks <- c(1,2,3,4)
lengths <- c(100,500)
alpha_val <- 0.05
nSims <- 1000

for(k in 1:length(Ks)){
  cat(paste0('\n - K ',Ks[k],' - '))
  results <- data.frame('k'=rep(Ks[k],length(lengths)),
                        'length'=lengths,
                        'TnSim'=NA, 'MnSim'=NA,'Mean'=NA)
  for(l in 1:length(lengths)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- readRDS(paste0(path,'MeanSin/',
                           'AltCase_MeanSinChange4_',Ks[k],
                           '_len',lengths[l],'_data.rds'))
    tmp_results <- data.frame('TnSim'=rep(NA,nSims), 'MnSim'=NA,'Mean'=NA)
    for(s in 1:nSims){
      cat(paste0(s,', '))

      set.seed(123*k+123^2*l+s)

      # Char
      TnMn <- detect_changepoint_final_TnAndMn(data[[s]], h=0,space='BM')
      # Tn
      Tn_sim <- TnMn[['Tn']]
      # Mn
      Mn_sim <- TnMn[['Mn']]

      mean_val <- mean_change(data[[s]],inc.pval = T)

      tmp_results[s,] <- c(Tn_sim$pval,
                           Mn_sim$pval,
                           mean_val[2])
    }

    saveRDS(tmp_results,paste0(path,'Mean/tmp/',
                               'mean_tmp_results_c',c,'_l',l,'.rds'))

  results[l,] <- c(changes[c],lengths[l],
                   colMeans(tmp_results<=alpha_val))
  saveRDS(results,paste0(path,'Mean/tmp/',
                         'mean_results_c',c,'.rds'))
  }

  saveRDS(results,paste0(path,'Mean/',
                         'mean_complete_c',c,'.rds'))
}
