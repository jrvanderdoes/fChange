####################################
#   KL - Loc Sim Change (50)
####################################
library(devtools); load_all()

change <- 0.15
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100)#,250,500)
nSims <- 1000
dat_res <- 50
sinKs <- c(1,5,10,25,50)

c <- 1
# for(c in 1:length(locs)){
for(K in 1:length(sinKs)){
  cat(paste0('\n - Loc ',locs[c],' - '))
  for(l in 1:length(lengths)){
    cat(paste0('\nLength ',lengths[l],': '))
    data <- list()
    for(s in 1:nSims){
      cat(paste0(s,','))
      set.seed(123*c+123^2*l+s+1)

      x_vals <- seq(0,1,length.out=dat_res)
      change <- sin(2*pi*sinKs[K]*x_vals)
      # Tn
      data[[s]] <- generate_data_fd(
        ns = c(lengths[l]*locs[c], lengths[l]-(lengths[l]*locs[c])),
        eigsList = list(c(1/sqrt(1:5)),
                        c(1/sqrt(1:5))),
        basesList = list(
          fda::create.bspline.basis(nbasis = 5, norder = 4)
        ),
        meansList = list(c(0),c(change)),
        distsArray = c("Normal",'Normal'),
        evals = seq(0, 1, length.out=dat_res),
        kappasArray = c(0),
        silent = T
      )
    }

    saveRDS(data,paste0(path,'KL/Loc/',
                        'AltCase_LocChange_sin_',sinKs[K],
                        '_len',lengths[l],'_data.rds'))
  }
}

####################################
#   KL - Loc Sim Change (50)
####################################
library(devtools); load_all()
path <- getwd()
path <- paste0(path,'/paperFiles/Alternative/KL/')

change <- 0.15
locs <- c(0.5,0.75,0.8,0.9,0.95)
lengths <- c(100)#,250,500)
nSims <- 1000
dat_res <- 50
sinKs <- c(1,5,10,25,50)

data <- data.frame('sinK'=sinKs,'locs'=NA,
                   'length'=NA, 'Tn'=NA,
                   'Mn'=NA, 'MeanMn'=NA,
                   'MeanMnResid'=NA, 'MeanTn'=NA,
                   'MeanTnResid'=NA)
c <- 1
for(K in 1:length(sinKs)){
  data_tmp <- readRDS(paste0(path,'Loc/tmp/sin/loc_sin_tmp_results_c',
                             c,K,'.rds'))
  data[K,-1]<- data_tmp
}

data_plot <- data %>% pivot_longer(cols = Tn:MeanTnResid)
ggplot() +
  geom_line(aes(x=sinK,y=value, group=name,color=name),
            data=data_plot) +
  theme_bw()
