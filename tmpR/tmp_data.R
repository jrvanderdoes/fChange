plot_fd(temperature)


tmp <- linear_imputatation(temperature$Sydney)
tmp <- na.omit(tmp)

plot_fd(tmp,FDReps = colnames(tmp))


msft_cidr <- readRDS('C:/Users/jerem/OneDrive/Documents/msft_cidr_30.rds')
