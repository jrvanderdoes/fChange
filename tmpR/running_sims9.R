#####################
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

  data <- readRDS(paste0('C:/Users/jerem/Downloads/fChangeData/',
                         'OUNull_len',l_idx,'_data.rds'))

  for(i in 1:nTrials){
    cat(paste0(i,', '))

    set.seed(123456*l_idx+i)
    tn0 <- detect_changepoint_final_Tn(data[[i]], h=0,space='BM')
    tn13 <- detect_changepoint_final_Tn(data[[i]], h=lengths[l_idx]^(1/3),space='BM')
    tnd13 <- detect_changepoint_final_Tn(data[[i]], h=2*lengths[l_idx]^(1/3),space='BM')

    results[i,] <- c(lengths[l_idx],
                     tn0$pval, tn13$pval, tnd13$pval)

    if((nTrials %% 100) == 0)
      saveRDS(results,paste0('C:/Users/jerem/Downloads/fChangeData/',
                             'Tn_OUNull_len',l_idx,'.rds'))
  }
}


#####################
## Mild Dependence
library(devtools);load_all()

# Normal
nTrials <- 1000
lengths <- c(100,250,500)
dat_res <- 21

for(l_idx in 1:length(lengths)){
  cat(paste0('\n',l_idx,': '))
  results <- data.frame('length'=rep(NA, nTrials),
                        'tnh0'=NA, 'tnh13'=NA, 'tnhd13'=NA)

  data <- readRDS(paste0('C:/Users/jerem/Downloads/fChangeData/',
                         'OUNull_len',l_idx,'_data_mild.rds'))

  for(i in 1:nTrials){
    cat(paste0(i,', '))

    set.seed(123456*l_idx+i)
    tn0 <- detect_changepoint_final_Tn(data[[i]], h=0,space='BM')
    tn13 <- detect_changepoint_final_Tn(data[[i]], h=lengths[l_idx]^(1/3),space='BM')
    tnd13 <- detect_changepoint_final_Tn(data[[i]], h=2*lengths[l_idx]^(1/3),space='BM')

    results[i,] <- c(lengths[l_idx],
                     tn0$pval, tn13$pval, tnd13$pval)

    if((nTrials %% 100) == 0)
      saveRDS(results,paste0('C:/Users/jerem/Downloads/fChangeData/',
                             'Tn_OUNull_len',l_idx,'_mild.rds'))
  }
}


#####################
## Strong Dependence
library(devtools);load_all()

# Normal
nTrials <- 1000
lengths <- c(100,250,500)
dat_res <- 21

for(l_idx in 1:length(lengths)){
  cat(paste0('\n',l_idx,': '))
  results <- data.frame('length'=rep(NA, nTrials),
                        'tnh0'=NA, 'tnh13'=NA, 'tnhd13'=NA)

  data <- readRDS(paste0('C:/Users/jerem/Downloads/fChangeData/',
                         'OUNull_len',l_idx,'_data_strong.rds'))

  for(i in 1:nTrials){
    cat(paste0(i,', '))

    set.seed(123456*l_idx+i)
    tn0 <- detect_changepoint_final_Tn(data[[i]], h=0,space='BM')
    tn13 <- detect_changepoint_final_Tn(data[[i]], h=lengths[l_idx]^(1/3),space='BM')
    tnd13 <- detect_changepoint_final_Tn(data[[i]], h=2*lengths[l_idx]^(1/3),space='BM')

    results[i,] <- c(lengths[l_idx],
                     tn0$pval, tn13$pval, tnd13$pval)

    if((nTrials %% 100) == 0)
      saveRDS(results,paste0('C:/Users/jerem/Downloads/fChangeData/',
                             'Tn_OUNull_len',l_idx,'_strong.rds'))
  }
}
