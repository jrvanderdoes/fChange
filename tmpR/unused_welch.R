path <- getwd()
path <- paste0(path,'/paperFiles/Null/KL/')
################################################
##  Welch - Tn 0.01
################################################
alp=0.01
###################
## No Dependence
library(devtools);load_all()

# Normal
nTrials <- 1000
lengths <- c(100,250,500)
dat_res <- 21

for(l_idx in 1:length(lengths)){
  cat(paste0('\n',l_idx,': '))
  results <- data.frame('length'=rep(NA, nTrials),
                        'tnh0'=NA, 'tnh13'=NA, 'tnhd13'=NA)

  data <- readRDS(paste0(path,'Data/',
                         'NullCase_len',l_idx,'_data.rds'))

  for(i in 1:nTrials){
    cat(paste0(i,', '))

    set.seed(123456*l_idx+i)
    tn0 <- detect_changepoint_Welch(
      X=data[[i]],
      h=0,space='BM',alpha=alp)
    tn13 <- detect_changepoint_Welch(
      X=data[[i]],
      h=lengths[l_idx]^(1/3),space='BM',alpha=alp)
    tnd13 <- detect_changepoint_Welch(
      X=data[[i]],
      h=2*lengths[l_idx]^(1/3),space='BM',alpha=alp)

    results[i,] <- c(lengths[l_idx],
                     tn0$detected, tn13$detected, tnd13$detected)

    if((nTrials %% 100) == 0)
      saveRDS(results,paste0(path,'tmpResults/', 'Tn_NullCase_welch01_len',l_idx,'.rds'))
  }

  saveRDS(results,paste0(path,'Results/', 'Tn_NullCase_welch01_len',l_idx,'.rds'))
}
