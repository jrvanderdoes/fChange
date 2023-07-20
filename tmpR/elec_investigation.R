library(tidyverse)
library(devtools)
library(fda)
load_all()

path <- 'C:/Users/jerem/Downloads/ElecData/Portugal/Data/'

elec_port <-
  do.call(rbind,
          lapply(paste0(path,dir(path)),read.csv,
                 sep = ';',skip=1,header = F))[,-c(6,7)]
colnames(elec_port) <- c('year','month','day','hour','price')
elec_port <- na.omit(elec_port)

elec_port <- elec_port %>%
  pivot_wider(id_cols = hour,
              names_from = c(year, month, day),
              values_from = price) %>%
  slice(-25) %>%
  as.data.frame() %>%
  linear_imputatation() %>%
  as.data.frame()

electricity_fd <- fda::Data2fd(
  argvals = 1:24,
  y = as.matrix(elec_port[,-1]),
  basisobj = create.bspline.basis(1:24))

# Play with more components
nPCs <- 10
electricity_fpca <- fda::pca.fd(electricity_fd,nharm = nPCs)
electricity_fpca_comp <- electricity_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for(i in 1:nPCs){
  ts_dat[[i]] <- ts(electricity_fpca_comp[,i], freq=7)
  #comps[[i]] <- forecast::ets(ts_dat[[i]])
  #comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

electricity_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
#   Want: 24 x 365
#     forecast: 365 x 3
#    coefs: 26  x 3
#      coefs %*% comps: 26 x 365
#    Eval: 24 x 26
#       eval %*% orig: 24 x 365
orig_coefs <- electricity_fpca$harmonics$coefs %*%
  t(electricity_fpca_forecast)
eval_fd_vals <- eval.basis(1:24,electricity_fd$basis) %*%
  orig_coefs


change_sc <- detect_changepoint_singleCov(
  X=eval_fd_vals, space = 'PC', h = 0,
  x = seq(0,1,length.out=50))

curveInt <- function(X){
  val <- rep(NA,ncol(X))
  for(i in 1:ncol(X)){
    val[i] <- pracma::trapz(seq(0,1,length.out=24),X[,i])
  }
  val
}

vals <- curveInt(eval_fd_vals)
vals_hm <- huskovaMeintanis(vals)

change_sc <- detect_changepoint_singleCov(
  X=eval_fd_vals, space = 'PC', h = 0,
  x = seq(0,1,length.out=50))

#####################################
eval_fd_vals24 <- eval.basis(1:24,electricity_fd$basis) %*%
  orig_coefs
eval_fd_vals2 <- eval.basis(seq(1,24,length.out=2),
                            electricity_fd$basis) %*%
  orig_coefs
eval_fd_vals100 <- eval.basis(seq(1,24,length.out=100),
                            electricity_fd$basis) %*%
  orig_coefs

compute_Tn(eval_fd_vals24,M = 20000)
compute_Tn(eval_fd_vals2,M = 20000)
compute_Tn(eval_fd_vals100,M = 20000)

W_tmp <- computeSpaceMeasuringVectors(20000,"RN",eval_fd_vals100)
compute_Tn(eval_fd_vals24,M = 5000,
           W = W_tmp[1:24,] )
compute_Tn(eval_fd_vals2,M = 5000,
           W = W_tmp[1:2,])
compute_Tn(eval_fd_vals100,M = 5000,
           W = W_tmp)

change_sc2 <- detect_changepoint_singleCov(
  X=eval_fd_vals2, space = 'RN', h = 0,
  x = seq(0,1,length.out=50)) # Change
change_sc24 <- detect_changepoint_singleCov(
  X=eval_fd_vals24, space = 'RN', h = 0,
  x = seq(0,1,length.out=50)) # No Change
