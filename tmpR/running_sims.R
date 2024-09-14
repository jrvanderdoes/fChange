library(devtools); load_all()

path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')
changes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
lengths <- c(100)#,250,500)
alpha_val <- 0.05
nSims <- 1000

for(c in 1:length(changes)){
  cat(paste0('\n - Change ',changes[c],' - '))
  results <- data.frame('change'=rep(changes[c],length(lengths)),
                        'length'=lengths,
                        'Mean'=NA,
                        'MeanResid'=NA,
                        'MeanTn'=NA,
                        'MeanTnResid'=NA)
  for(l in 1:length(lengths)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- readRDS(paste0(path,'Mean/',
                           'AltCase_MeanChange_',changes[c],
                           '_len',lengths[l],'_data.rds'))
    tmp_results <- data.frame('Mean'=rep(NA,nSims),
                              'MeanResid'=NA,
                              'MeanTn'=NA,
                              'MeanTnResid'=NA)
    for(s in 1:nSims){
      cat(paste0(s,', '))
      tmp_dat <- data[[s]] -
        data.frame(matrix(0,nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2),
                   matrix(changes[c],nrow=nrow(data[[s]]),
                          ncol=ncol(data[[s]])/2))

      set.seed(123*c+123^2*l+s)
      tmp <- mean_change(data[[s]], h=0)
      tmp1 <- mean_change(data[[s]],data1=tmp_dat, h=0)
      tmpTn <- mean_change_Tn(data[[s]], h=0)
      tmpTn1 <- mean_change_Tn(data[[s]],data1=tmp_dat, h=0)
      tmp_results[s,] <- c(tmp,tmp1,tmpTn,tmpTn1)
    }

    saveRDS(tmp_results,paste0(path,'Mean/newTry/',
                               'dist_tmp_results_c',c,'_l',l,'_meanMethod.rds'))

    results[l,] <- c(changes[c],lengths[l],colMeans(!is.na(tmp_results)))
    saveRDS(results,paste0(path,'Mean/newTry/',
                           'mean_results_c',c,'_meanMethod.rds'))
  }

  saveRDS(results,paste0(path,'Mean/newTry/Final/',
                         'dist_complete_c',c,'_meanMethod.rds'))
}
