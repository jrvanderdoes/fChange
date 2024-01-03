set.seed(123456)
res <- 20
M <- 50

hVal <- 28

.generateFunctionalIID <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = FALSE)

  t(fiid)
}

x <- seq(0, 1, length.out = res)
MJ <- M * length(x)

obs_opts <- c(100, 500, 1000, 2500, 5000)

result <- data.frame('pval100'=NA, 'pval500'=NA,
                     'pval1000'=NA, 'pval2500'=NA,
                     'pval5000'=NA)
nSims <- 200

set.seed(12345)
for(j in 1:nSims){
  cat(paste0(j,', '))
  data_use <- .generateFunctionalIID(res,max(obs_opts))

  for(i in 1:length(obs_opts)){

    sqrtMat <- .compute_sqrtMat(Cov_M=M, space='BM', X=data_use[,1:obs_opts[i]],
                                x=x,h=hVal,K=bartlett_kernel)

    Cov_iid <- detect_changepoint_singleCov(
      X=data_use[,1:500],x = x,TN_M = M,K = bartlett_kernel,
      nSims = 1000, h=hVal, space='BM', Cov_M=M, silent = T,
      sqrtMat = sqrtMat)

    result[j,i] <- Cov_iid$pval
  }
  #saveRDS(result,'C:/Users/jerem/Downloads/obs_run_sqrtMat_updateCovMat.rds')
  saveRDS(result,'C:/Users/jerem/Downloads/obs_run_sqrtMat_updatef.rds')
}

result_plot <- result %>% pivot_longer(cols=pval100:pval5000)
result_plot$name<-factor(result_plot$name,levels=paste0('pval',obs_opts))
levels(result_plot$name) <- c('100', '500', '1000', '2500', '5000')

ggplot(result_plot,aes(x=name,y=value)) +
  geom_boxplot() +
  geom_jitter(alpha=0.25) +
  geom_hline(aes(yintercept=0.05)) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  ylab('p-value') +
  xlab('Cov Estimate Length') +
  geom_text(data=data.frame(),
            aes(x=c('100', '500', '1000', '2500', '5000'),
                y=1.1,
                label=as.numeric(colMeans(result<=0.05))),
            size=4) +
  ggtitle('500 Length TS (200 Sims)')

