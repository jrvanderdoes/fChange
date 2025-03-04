rates <- read.csv('C:/Users/jerem/Downloads/yield-curve-rates-1990-2021.csv')
rates <- t(rates)
months <- c(1,2,3,4,6,12,24,36,5*12,7*12,10*12,20*12,30*12)

rates <- cbind(data.frame(c('Date',months)),
               rates)
rates_only <- dplyr::mutate_all(rates[-1,-1], as.numeric)
rates_only <- linear_imputatation(data = rates_only,
                                  evalPts = months)

#changes <- detect_changepoint_singleCov(as.matrix(house[,-c(1,ncol(house))]))
mn <- compute_Mn(X = rates_only)
#changes <-


plot_fd(data = rates_only,
        changes=mn$location,
        curve_points = months)
plot_fd(data = rates_only,
        curve_points = months,
        FDReps = as.Date(as.character(as.vector(rates[1,])),'%m/%d/%y'))

compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}

rates_cidr <- compute_cidr(rates_only)
mn_cidr <- compute_Mn(X = rates_cidr)

plot_fd(data = rates_cidr,changes=mn_cidr$location,
        curve_points = months)
plot_fd(data = rates_cidr,
        curve_points = months,
        FDReps = as.Date(as.character(as.vector(rates[1,])),'%m/%d/%y'))
