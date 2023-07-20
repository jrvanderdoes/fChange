
data_length <- 50
nSims <- 200
dep <- 0

result <- data.frame("Cov_iid"=rep(NA,nSims),
                     "Cov_fn5"=NA,
                     "Cov_f2n5"=NA,
                     "Perm_iid"=NA,
                     "Perm_fn5"=NA,
                     "Perm_f2n5"=NA,
                     "Welch_iid"=NA,
                     "Welch_fn5"=NA,
                     "Welch_f2n5"=NA)
for(i in 1:nSims){
  cat(paste0(i,', '))
  null_case <- generate_data_fd(ns = c(data_length),
                                eigsList = list(c(3,2,1,0.5)),
                                basesList = list(fda::create.bspline.basis(nbasis=4, norder=4)),
                                meansList = c(0),
                                distsArray = c('Normal'),
                                evals = seq(0,1,0.05),
                                kappasArray = c(dep),
                                silent = T)
  Cov_iid <- detect_changepoint_singleCov(
    X=null_case, nSims = 1000, h=0,
    space='BM', Cov_M=50, silent = T)
  Cov_n5 <- detect_changepoint_singleCov(
    X=null_case, nSims = 1000, h=floor(data_length^(1/5)),
    space='BM', Cov_M=50, silent = T)
  Cov_2n5 <- detect_changepoint_singleCov(
    X=null_case, nSims = 1000, h=floor(2*data_length^(1/5)),
    space='BM', Cov_M=50, silent = T)

  Perm_iid <- generalized_resampling(
    X=null_case, blockSize = 1, fn = compute_Tn,
    space='BM', M=1000,silent = T)
  Perm_n5 <- generalized_resampling(
    X=null_case, blockSize = floor(data_length^(1/5)),
    fn = compute_Tn, space='BM', M=1000,silent = T)
  Perm_2n5 <- generalized_resampling(
    X=null_case, blockSize = floor(2*data_length^(1/5)),
    fn = compute_Tn, space='BM', M=1000,silent = T)

  Approx_iid <- welch_approximation(
    X=null_case, alpha = 0.05, TVal = ncol(null_case),
    W = NULL, W1 = NULL, M=1000, h = 0,
    K = bartlett_kernel)
  Approx_n5 <- welch_approximation(
    X=null_case, alpha = 0.05, TVal = ncol(null_case),
    W = NULL, W1 = NULL, M=1000, h = floor(data_length^(1/5)),
    K = bartlett_kernel)
  Approx_2n5 <- welch_approximation(
    X=null_case, alpha = 0.05, TVal = ncol(null_case),
    W = NULL, W1 = NULL, M=1000, h = floor(2*data_length^(1/5)),
    K = bartlett_kernel)

  Tn <- mean(Cov_iid$value,Cov_n5$value,Cov_2n5$value)

  result[i, ] <- c(Cov_iid$pval <= 0.05,
                   Cov_n5$pval <= 0.05,
                   Cov_2n5$pval <= 0.05,

                   Perm_iid$pval <= 0.05,
                   Perm_n5$pval <= 0.05,
                   Perm_2n5$pval <= 0.05,

                   Approx_iid <= Tn,
                   Approx_n5 <= Tn,
                   Approx_2n5 <= Tn)
}

colSums(result,na.rm = T)/nSims
