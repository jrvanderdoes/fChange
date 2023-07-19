
library(fda)


electricity <- fChange::electricity
electricity_fd <- fda::Data2fd(1:24,electricity)
# sum(fda::eval.fd(1:24,electricity_fd)-electricity)
#   eval.fd: eval.basis(1:24,electricity_fd$basis) %*% electricity_fd$coefs

# Play with more components
nPCs <- 5
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
orig_coefs <- electricity_fpca$harmonics$coefs %*% t(electricity_fpca_forecast)
eval_fd_vals <- eval.basis(1:24,electricity_fd$basis) %*% orig_coefs

# Check Status
#tmp_acf <- fdaACF::obtain_FACF(t(eval_fd_vals),v=seq(0,1,length.out=24),nlags=20)

# CE Method (1-dim)
tmp_dat <- apply(eval_fd_vals,MARGIN = 2,
      function(x,l){rep(mean(x),l)},
      l=nrow(eval_fd_vals))
ce_est <- ecp::e.cp3o(t(tmp_dat))
plot(tmp_dat[1,])
abline(v = ce_est$estimates)

tmp_dat_c <- detect_changepoint_singleCov(
                    X=tmp_dat, space = 'RN')
tmp_dat_c$pval
tmp_dat_c1 <- detect_changepoint_singleCov(
                    X=eval_fd_vals, space = 'RN')
tmp_dat_c1$pval
# plot_fd(tmp_dat, compute_Mn(tmp_dat, which.Mn=TRUE))


# True Data
tmp_mc <- mean_change(eval_fd_vals,inc.pval = T)
tmp_mc
tmp_cc <- cov_change(eval_fd_vals)
tmp_cc

tmp_sc <- detect_changepoint_singleCov(X=eval_fd_vals, h=0)
tmp_sc$pval
compute_Mn(eval_fd_vals, which.Mn=TRUE)

tmp_gs <- generalized_resampling(X=eval_fd_vals, blockSize=1,
                                 fn=compute_Tn, iters=1000,
                                 replace=FALSE)
tmp_gs

tmp_wa <- welch_approximation(eval_fd_vals, h=0)
tmp_wa

## Try Huskova Paper code with integrate data


