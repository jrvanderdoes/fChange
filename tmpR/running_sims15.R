path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/FAR1/')

####################################
#   FAR(1) - Mean Change
####################################
library(devtools); load_all()

changes <- c(0, 0.15, 0.25, 0.3, 0.35, 0.4)
lengths <- c(100,250,500)
nSims <- 1000
dat_res <- 50
alpha_val <- 0.05

c <- 4
#for(c in 1:length(changes)){
results <- data.frame('changes'=rep(changes[c],length(lengths)),
                      'length'=lengths,
                      'Tn'=NA,'Mn'=NA,'Mean'=NA)

l <- 3

cat(paste0('\nLength ',lengths[l],': '))
data <- readRDS(
  paste0(path,'Mean1/',
         'AltCase_MeanChange_',changes[c],
         '_len',lengths[l],'_FAR1data.rds'))
tmp_results <- data.frame('TnSim'=rep(NA,nSims),
                          'MnSim'=NA,'Mean'=NA)
for(s in 701:800){
  cat(paste0(s,', '))

  set.seed(123*c+123^2*l+s)

  # Chars
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

# saveRDS(tmp_results,
#         paste0(path,'Mean1/tmp/tmp/dat_tmp.rds',changes[c],
#                '_len',lengths[l],'_FAR1data.rds'))
