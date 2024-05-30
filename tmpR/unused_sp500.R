library(tidyverse)
library(ggplot2)

##########################
path <- "C:/Users/jerem/Downloads/newstock/"
stock <- 'SP500'
years <- 2019:2022#2019:2023
data<- list()

for(i in 1:length(years)){
  data[[2*(i-1)+1]] <- read.csv(paste0(path,stock,years[i],"M.csv"))
  data[[2*(i-1)+2]] <- read.csv(paste0(path,stock,years[i],"M1.csv"))
}
data <- do.call("rbind",data)
data$Date <- as.POSIXct(data$Local.time,format='%d.%m.%Y %H:%M:%S')
data$YearMonth <- paste0(year(data$Date),'_',month(data$Date))
data$YearMonthDay <- paste0(year(data$Date),'_',month(data$Date),'_',day(data$Date))
data$Day <- day(data$Date)
data$Hour <- hour(data$Date)
data$Min <- minute(data$Date)

data_Daily <- data.frame(matrix(NA,
                                ncol=length(unique(data$YearMonthDay)),
                                nrow=24*60))
colnames(data_Daily) <- unique(data$YearMonthDay)
times <- data[data$YearMonthDay==colnames(data_Daily)[1], "Date"]
times <- paste0(hour(times),'_',minute(times),'_',second(times))
for(i in 1:ncol(data_Daily)){
  dat0 <- data.frame(
    'Date'=as.POSIXct(paste(colnames(data_Daily)[i], times),
                      format = "%Y_%m_%d %H_%M_%S"))
  obs_data <- data[data$YearMonthDay==colnames(data_Daily)[i],c('Date',"Close")]

  dat <- merge(dat0,obs_data, by='Date',
               all.x = T)
  if(nrow(dat)==24*60){
    data_Daily[,i] <- dat$Close
  } else{
    dat <- merge(dat0, aggregate(Close ~ Date, dat, mean), by='Date',
                 all.x = T)
    data_Daily[,i] <- dat$Close
  }
}

# Clean
#data_Daily <- data_Daily[(9.5*60):(16*60),]  # Remove non-trading hrs
data_Daily <- data_Daily[(10*60):(16*60),]  # Remove non-trading hrs
uniquelength <- sapply(data_Daily,function(x) length(unique(x)))
data_Daily <- subset(data_Daily, select=uniquelength>1) # Remove non-trading days
dates <- colnames(data_Daily)
data_Daily <- linear_imputatation(data_Daily,use.prev.curve = T)
colnames(data_Daily) <- dates

# plot_fd(data_Daily)


# bs <- binary_segmentation(X=data_Daily,
#                           statistic = 'Tn',
#                           method = 'Sim') # all
# plot_fd(data_Daily,bs)


data_cidr <- data_Daily
for(i in 1:nrow(data_Daily)){
  data_cidr[i,] <- 100*( log(data_Daily[i,]) - log(data_Daily[1, ]))
}
#plot_fd(data_cidr)

bs_cidr <- binary_segmentation(data_cidr,'Tn','Sim')

plot_fd(data_cidr,bs_cidr)
