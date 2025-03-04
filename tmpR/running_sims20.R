path <- getwd()
path <- paste0('C:\\Users\\jerem\\Downloads\\PlayHere\\')

####################################
#   Mean Change - CE Method
####################################
library(devtools); load_all()

# changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
Ks <- c(1,2,3,4)
lengths <- c(100,500)
alpha_val <- 0.05
nSims <- 1000

k <- 4
#for(k in 1:length(Ks)){
cat(paste0('\n - K ',Ks[k],' - '))
results <- data.frame('k'=rep(Ks[k],length(lengths)),
                      'length'=lengths,
                      'TnSim'=NA, 'MnSim'=NA,'Mean'=NA)
l = 1
# for(l in 1:1){
cat(paste0('\nLength ',lengths[l],': '))
data <- readRDS(paste0(path,'MeanSin/',
                       'AltCase_MeanSinChange4_',Ks[k],
                       '_len',lengths[l],'_data.rds'))
tmp_results <- data.frame('TnSim'=rep(NA,nSims), 'MnSim'=NA,'Mean'=NA)
for(s in 500:1){
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

# saveRDS(tmp_results,paste0(path,'Mean/tmp/',
#                            'mean_tmp_results_c',c,'_l',l,'.rds'))

# results[l,] <- c(changes[c],lengths[l],
#                  colMeans(tmp_results<=alpha_val))
# saveRDS(results,paste0(path,'Mean/tmp/',
#                        'mean_results_c',c,'.rds'))
# }

#   saveRDS(results,paste0(path,'Mean/',
#                          'mean_complete_c',c,'.rds'))
# }


##################################
# saveRDS(tmp_results,
#         paste0(path,'MeanSin/tmp/dat_tmp',Ks[k],
#                '_len',lengths[l],'.rds'))



###############################
# tmp_results_old <- tmp_results
# tmp1 <- readRDS(paste0(path,'MeanSin/tmp/dat_tmp',Ks[k],
#                  '_len',lengths[l],'.rds'))
# idxs <- 501:1000
# tmp_results[idxs,] <- tmp1[idxs,]
# colMeans(tmp_results)
#
# saveRDS(tmp_results,
#         paste0(path,'MeanSin/tmp/tmp_results',Ks[k],
#                        '_len',lengths[l],'.rds'))

###############################
# results <- data.frame(expand.grid('K'=Ks, 'length'=c(100,500)),
#                       'Tn'=NA,'Mn'=NA,'Mean'=NA)
# #results <- readRDS(paste0(path,'MeanSin/complete_results10.rds'))
# results[results$K==Ks[k] &
#           results$length==lengths[l], ] <-
#   c(Ks[k],lengths[l],colMeans(tmp_results<=0.05))
# results
# saveRDS(results,paste0(path,'MeanSin/complete_results2.rds'))

tmp_results',Ks[k], '_len',lengths[l],'.rds'
