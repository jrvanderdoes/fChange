## ELECTRICITY BS
library(fda)

electricity <- fChange::electricity
electricity_fd <- fda::Data2fd(1:24,electricity)

nPCs <- 3
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

#
# generalized_resampling(X=eval_fd_vals,
#                        blockSize=1,
#                        fn=compute_Mn, iters=1000, replace=F)


changes=c(94,202,275)
pdf("electricity_resids.pdf")
plot_fd(eval_fd_vals,changes = changes,interactive = FALSE,
        val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()

pdf("electricity.pdf")
plot_fd(electricity,changes = changes,interactive = FALSE,
        val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()

