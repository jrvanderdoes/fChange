path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/OU/')

####################################
#   Mean Change
####################################
library(devtools); load_all()

changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4)
lengths <- c(100,250,500)
alpha_val <- 0.05
nSims <- 1000

# for(c in 1:length(changes)){

c <- 8

cat(paste0('\n - Change ',changes[c],' - '))
results <- data.frame('change'=rep(changes[c],length(lengths)),
                      'length'=lengths,
                      'TnSim'=NA,
                      'MnSim'=NA,'Mean'=NA)
l <- 3
# for(l in 1:length(lengths)){
cat(paste0('\nLength ',lengths[l],': '))
data <- readRDS(paste0(path,'Mean1/',
                       'AltCase_MeanChange_',changes[c],
                       '_len',lengths[l],'_OUdata.rds'))
tmp_results <- data.frame('TnSim'=rep(NA,nSims),
                          'MnSim'=NA,'Mean'=NA)
for(s in 551:700){
  cat(paste0(s,', '))

  set.seed(123*c+123^2*l+s)

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

#     saveRDS(tmp_results,paste0(path,'Mean/tmp/',
#                                'meanOU_tmp_results_c',c,'_l',l,'.rds'))
#
#     results[l,] <- c(changes[c],lengths[l],colMeans(tmp_results<=alpha_val))
#     saveRDS(results,paste0(path,'Mean/tmp/',
#                            'meanOU_results_c',c,'.rds'))
#   }
#
#   saveRDS(results,paste0(path,'Mean/', 'meanOU_complete_c',c,'.rds'))
# }
#
# saveRDS(tmp_results,
#         paste0(path,'Mean1/tmp/tmp/dat_tmp.rds',changes[c],
#        '_len',lengths[l],'_FAR1data.rds'))
