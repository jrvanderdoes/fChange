ress <- c(20,40,60,80)
Ms <- c(20,40,60,80)
nSims <- 250

hVal <- 0

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

result <- data.frame('sim'=rep(NA,nSims*length(ress)*length(Ms)),'res'=NA,'M'=NA,'pval'=NA)

idx <- 0
for(i in 1:nSims){
  data_use <- .generateFunctionalIID(resolution = 20,N = 500)

  for (res in ress){
    for (M in Ms){
      idx <- idx + 1
      cat(paste0(idx, ', '))
      Cov_iid <- detect_changepoint_singleCov_n(
        X=data_use,x = seq(0, 1, length.out = res),
        TN_M = M,K = bartlett_kernel,
        nSims = 1000, h=hVal, space='BM', Cov_M=M, silent = T)

      result[idx,] <- c(i,res,M,Cov_iid$pval)
    }
  }
}
#saveRDS(result,'C:\\Users\\jerem\\Downloads\\ShowFigures\\img1.rds')
result_plot <- result
result_plot$resM <- paste0(result_plot$res,result_plot$M)
result_plot$label <- NA
for (res in ress){
  for (M in Ms){
    idx <- idx+1
    result_plot[result_plot$res==res &
                  result_plot$M == M,'label'] <-
      c(mean(result_plot[result_plot$res==res & result_plot$M==M,'pval']<=0.05),
        rep(NA,nrow(result_plot[result_plot$res==res &
                                  result_plot$M == M,])-1))
  }
}

ggplot(result_plot) +
  facet_grid(rows=vars(res),cols=vars(M)) +
  geom_boxplot(aes(y=pval)) +
  geom_jitter(aes(x=0,y=pval),alpha=0.25) +
  geom_hline(aes(yintercept=0.05)) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_text(aes(x=0,
                y=1.1,
                label=label),
            size=4) +
  ylab('p-value (res)') +
  xlab('(M)') +
  theme(axis.text.x = element_blank()) +
  ggtitle('500L TS with res 20 (250sim,h=0,iid, nmatch) - Cov Matrix')
