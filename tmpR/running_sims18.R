####################################
#   Eigen Change - CE Method
####################################
library(devtools); load_all()

path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')

changes <- c(1, 2, 4, 6, 8)
lengths <- c(100,500)
alpha_val <- 0.05
nSims <- 1000

c <- 5
cat(paste0('\n - Change ',changes[c],' - '))
results <- data.frame('change'=rep(changes[c],length(lengths)),
                      'length'=lengths,
                      'TnSim'=NA,
                      'MnSim'=NA,
                      'Cov'=NA)

l <- 2
cat(paste0('\nLength ',lengths[l],': '))
data <- readRDS(paste0(path,'Eigen/',
                       'AltCase_EigenChange_',changes[c],
                       '_len',lengths[l],'_data1.rds'))
tmp_results <- data.frame('TnSim'=rep(NA,nSims),
                          'MnSim'=NA,
                          'Cov'=NA)

sims <- 801:1000
for(s in sims){
  cat(paste0(s,', '))

  set.seed(123*c+123^2*l+s)
  # Char
  TnMn <- detect_changepoint_final_TnAndMn(data[[s]], h=0,space='BM')
  # Tn
  Tn_sim <- TnMn[['Tn']]
  # Mn
  Mn_sim <- TnMn[['Mn']]

  tryCatch({
    cov_sim <- cov_change(data[[s]],kappa=0)
  },
  error=function(e){
    # This would shift seed enough perhaps to solve issue
    cov_sim <- cov_change(data[[s]],kappa=0)
  })

  tmp_results[s,] <- c(Tn_sim$pval,
                       Mn_sim$pval,
                       cov_sim)
}

saveRDS(tmp_results,paste0(path,'/Eigen/500k0/c',
                           c,'_',5,'_500.rds'))

#     saveRDS(tmp_results,paste0(path,'Eigen/tmp1/',
#                                'tmp_results_c',c,'_l',l,'.rds'))
#
#     results[l,] <- c(changes[c],lengths[l],colMeans(tmp_results<=alpha_val))
#     saveRDS(results,paste0(path,'Eigen/tmp1/',
#                            'results_c',c,'.rds'))
#   }
#
#   saveRDS(results,paste0(path,'Eigen/tmp1/', 'complete_c',c,'.rds'))
# }
