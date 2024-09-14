path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/')

####################################
#   KL - Eigen Change (50)
####################################
library(devtools); load_all()

changes <- x(1,4,8,16)#c(1, 2, 4, 6, 8, 16)
lengths <- c(100,500)
nSims <- 1000
dat_res <- 50

for(c in 1:length(changes)){
  cat(paste0('\n - Change ',changes[c],' - '))
  for(l in 1:length(lengths)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- list()
    for(s in 1:nSims){
      cat(paste0(s,','))
      set.seed(123*c+123^2*l+s+1)
      # Tn
      data[[s]] <- generate_data_fd(
        ns = c(lengths[l]/2, lengths[l]/2),
        eigsList = list(c(1/sqrt(1:5)),
                        c(1/sqrt(1:5)*changes[c])),
        basesList = list(
          fda::create.bspline.basis(nbasis = 5, norder = 4)
        ),
        meansList = c(0,0),
        distsArray = c("Normal"),
        evals = seq(0, 1, length.out=dat_res),
        kappasArray = c(0),
        silent = T
      )
    }

    saveRDS(data,paste0(path,
                        'KL/Eigen/AltCase_EigenChange_',changes[c],
                        '_len',lengths[l],'_data2.rds'))
  }
}
