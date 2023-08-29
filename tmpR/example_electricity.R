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


elec_bs <- complete_binary_segmentation(
  data=eval_fd_vals,
  # test_statistic_function = compute_Tn,
  # changepoint_function = detect_changepoint_singleCov,
  trim_function=function(data, ...){max(2, floor(log(ncol(as.data.frame(data)))),na.rm=T)},
  changepoint_function=function(data, alpha, ...){
    func_val <- generalized_resampling(X=data,
        blockSize=1,
        fn=compute_Mn, iters=1000, replace=F, ...)

    ifelse(func_val$pval < alpha,
           compute_Mn(data, which.Mn=TRUE, ...),
           NA)
  },
  M=1000,
  final_verify = T,
  silent = F,
  space='BM', h=0)

# Paper Plotting
plot_fd(electricity, showticklabels = F,
        FD_axis_title = '',val_axis_title = '',
        eval_axis_title = '',
        eye = list(x = -.25, y = -1.5, z = 0.5),
        aspectratio=list(x=1.5,y=0.5,z=0.75))
plot_fd(eval_fd_vals,elec_bs, showticklabels = F,
        FD_axis_title = '',val_axis_title = '',
        eval_axis_title = '',
        eye = list(x = -.25, y = -1.5, z = 0.5),
        aspectratio=list(x=1.5,y=0.5,z=0.75))
plot_fd(electricity,elec_bs, showticklabels = F,
        FD_axis_title = '',val_axis_title = '',
        eval_axis_title = '',
        eye = list(x = -.25, y = -1.5, z = 0.5),
        aspectratio=list(x=1.5,y=0.5,z=0.75))

## Load saved
cps=c(94,202,275)
pdf("electricity_resids.pdf")
plot_fd(eval_fd_vals,CPs = cps,interactive = FALSE,
        val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()

pdf("electricity.pdf")
plot_fd(electricity,CPs = cps,interactive = FALSE,
        val_axis_title = '',res_axis_title = '',FD_axis_title = '',showticklabels = F)
dev.off()
