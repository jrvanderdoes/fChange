load_all()

nSims <- 1000
length <- 100
Ms <- c(5,10,20,30,35)
data <- list()

for(M in 1:length(Ms)){
  results_tmp <- data.frame('Tn'=rep(NA,nSims),'Tn95'=NA,
                            'Mn'=NA,'Mn95'=NA)
  for(s in nSims:1){
    cat(paste0(s,','))
    set.seed(12345*M+s+1)

    data[[s]] <- generate_data_fd(
      ns = c(length),
      eigsList = list(c(1/sqrt(1:5))),
      basesList = list(
        fda::create.bspline.basis(nbasis = 5, norder = 4)
      ),
      meansList = c(0),
      distsArray = c("Normal"),
      evals = seq(0, 1, length.out=50),
      kappasArray = c(0),
      silent = T
    )
    # Tn
    tmp <- detect_changepoint_final_TnAndMn(
      X = data[[s]], M=Ms[M], J=50, h=0)

    results_tmp[s,] <- c(tmp$Tn$value,quantile(tmp$Tn$gamProcess,probs = 0.95)[[1]],
                         tmp$Mn$value,quantile(tmp$Mn$gamProcess,probs = 0.95)[[1]])
  }

  saveRDS(results_tmp,
          paste0('C:\\Users\\jerem\\Downloads\\PlayHere\\M',
                 res, '_100.rds'))
}
