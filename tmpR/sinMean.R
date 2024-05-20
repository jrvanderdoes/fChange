####################################
#   Sin Mean
####################################
library(devtools); load_all()
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')


change <- 0.15
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100)
nSims <- 1000
dat_res <- 50
alpha_val <- 0.05
nSims <- 1000
dat_res <- 50
sinKs <- c(1,5,10,25,50)

c <- 1

K <- 1
# for(K in 1:length(sinKs)){
  cat(paste0('\n - Change ',sinKs[K],' - '))
  results <- data.frame('locs'=rep(sinKs[K],length(lengths)),
                        'length'=lengths,
                        'TnSim'=NA,'TnBoot'=NA,'TnWelch'=NA,
                        'MnSim'=NA,'MnBoot'=NA,'Mean'=NA)
  for(l in 1:length(lengths)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- readRDS(paste0(path,'Loc/',
                           'AltCase_LocChange_sin_',sinKs[K],
                           '_len',lengths[l],'_data.rds'))
    tmp_results <- data.frame('TnSim'=rep(NA,nSims),'TnBoot'=NA,
                              'TnWelch'=NA,
                              'MnSim'=NA,'MnBoot'=NA,'Mean'=NA)
    for(s in 1:nSims){
      cat(paste0(s,', '))

      set.seed(123*c+123^2*l+s)
      # Char
      TnMn <- detect_changepoint_final_TnAndMn(data[[s]], h=0,space='BM')
      # Tn
      Tn_sim <- TnMn[['Tn']]
      # Mn
      Mn_sim <- TnMn[['Mn']]

      # Mean
      tmp_dat <- data[[s]] -
        data.frame(matrix(0,nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2),
                   matrix(sin(2*pi*sinKs[K]*seq(0,1,length.out=50)),
                          nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2))
      mean_val <- mean_change(data[[s]],inc.pval = T)
      mean_val1 <- mean_change(data[[s]],data1=tmp_dat,inc.pval = T)
      mean_valTn <- mean_change_Tn(data[[s]],inc.pval = T)
      mean_valTn1 <- mean_change_Tn(data[[s]],data1=tmp_dat,inc.pval = T)

      mean_val <- mean_change(data[[s]],inc.pval = T)
      tmp_results[s,] <- c(Tn_sim$pval, Mn_sim$pval,
                           mean_val[2],
                           mean_val1[2],
                           mean_valTn[2],
                           mean_valTn1[2])
    }

    saveRDS(tmp_results,paste0(path,'Loc/tmp/sin/',
                               'loc_sin_tmp_results_c',c,'_l',l,'.rds'))

    results[l,] <- c(locs[c],lengths[l],colMeans(tmp_results<=alpha_val))
    saveRDS(results,paste0(path,'Loc/tmp/sin/',
                           'loc_sin_tmp_results_c',c,K,'.rds'))
  }

#   saveRDS(results,paste0(path,'Loc/', 'loc_sin_complete_c',c,'.rds'))
# }
