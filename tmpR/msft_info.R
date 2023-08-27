msft <- readRDS('C:/Users/jerem/OneDrive/Documents/msft_base.rds')
msft_use <- linear_imputatation(as.data.frame(msft[-1]))

msft_cidr <- msft_use
for(i in 1:nrow(msft_use)){
  msft_cidr[i,] <- 100*(log(msft_use[i,]) - log(msft_use[1,]))
}

plot_fd(msft_cidr)


#############

msft_cidr <- readRDS('C:/Users/jerem/OneDrive/Documents/msft_cidr.rds')

result_cidr_dc <- detect_changepoint_singleCov(msft_cidr, maxM=2500)
result_cidr_dc$pval

result_cidr_gr <- generalized_resampling(X=msft_cidr,
                                         blockSize = ncol(msft_cidr)^(1/3),
                                         fn=compute_Tn, iters=2500,
                                         replace = F)
result_cidr_gr
