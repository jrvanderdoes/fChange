####################################
#   FAR(1) - Mean Change (50)
####################################
library(devtools); load_all()
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/')


changes <- c(0, 0.15, 0.25, 0.3, 0.35, 0.4)
lengths <- c(100,250,500)
nSims <- 1000
dat_res <- 50

for(c in 1:length(changes)){
  cat(paste0('\n - Change ',changes[c],' - '))
  for(l in 1:length(lengths)){
    if(l==2) next
    cat(paste0('\nLength ',lengths[l],': '))
    data <- list()
    for(s in 1:nSims){
      cat(paste0(s,','))
      set.seed(123*c+123^2*l+s)
      # Tn
      data[[s]] <- generateFAR1(resolution = dat_res, d=0,
                                N = lengths[l],burnin = 1000)
      data[[s]][,(lengths[l]/2 + 1):lengths[l]] <-
        data[[s]][,(lengths[l]/2 + 1):lengths[l]] + changes[c]
    }

    saveRDS(data,paste0(path,'FAR1/Mean1/',
                        'AltCase_MeanChange_',changes[c],
                        '_len',lengths[l],'_FAR1data.rds'))
  }
}
