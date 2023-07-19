compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}

ford <- readRDS('C:/Users/jerem/Downloads/ford.rds')
ctn <- readRDS('C:/Users/jerem/Downloads/ctn.rds')
goog <- readRDS('C:/Users/jerem/Downloads/goog.rds')
msft <- readRDS('C:/Users/jerem/Downloads/msft.rds')

plot_fd(as.data.frame(ford[-1]))
# plot_fd(compute_cidr(linear_imputatation(as.data.frame(ford[-1]))))
plot_fd(linear_imputatation(as.data.frame(ctn[-1])))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(ctn[-1]))))
plot_fd(as.data.frame(goog[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(goog[-1]))))
plot_fd(as.data.frame(msft[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(msft[-1]))))

spx <- readRDS('C:/Users/jerem/Downloads/spx.rds')
tsla <- readRDS('C:/Users/jerem/Downloads/tsla.rds')
cake <- readRDS('C:/Users/jerem/Downloads/cake.rds')
nflx <- readRDS('C:/Users/jerem/Downloads/nflx.rds')
oil <- readRDS('C:/Users/jerem/Downloads/oil.rds')

plot_fd(as.data.frame(spx[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(spx[-1]))))
plot_fd(as.data.frame(tsla[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(tsla[-1]))))
plot_fd(as.data.frame(cake[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(cake[-1]))))
plot_fd(as.data.frame(nflx[-1]))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(nflx[-1]))))
plot_fd(linear_imputatation(as.data.frame(oil[-1])))
#plot_fd(compute_cidr(linear_imputatation(as.data.frame(oil[-1]))))


tmp1 <- detect_changepoint(linear_imputatation(as.data.frame(oil[-1])))
tmp1$pval

tmp2 <- generalized_resampling(X=linear_imputatation(as.data.frame(oil[-1])),
    blockSize=ncol(linear_imputatation(as.data.frame(oil[-1])))^(1/3),
    fn=compute_Tn, iters=1000, replace=FALSE)
tmp2
compute_Tn(linear_imputatation(as.data.frame(oil[-1])))

results <- .autoElbow_method(data=linear_imputatation(as.data.frame(oil[-1])))
plot_fd(data=linear_imputatation(as.data.frame(oil[-1])), CPs=results)
