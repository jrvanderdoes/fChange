##########################
BTC <- read.csv("C:/Users/jerem/Downloads/BTC-USD.csv")
BTC$Date <- as.Date(BTC$Date)
BTC$YearMonth <- paste0(substr(BTC$Date,1,4),'_',substr(BTC$Date,6,7))
BTC$Day <- substr(BTC$Date,9,10)

BTC_Monthly <- data.frame(matrix(NA,ncol=length(unique(BTC$YearMonth)),
                                 nrow=30))
colnames(BTC_Monthly) <- unique(BTC$YearMonth)
for(i in 1:ncol(BTC_Monthly)){
  obs_data <- BTC[BTC$YearMonth==colnames(BTC_Monthly)[i],"Adj.Close"]
  tmp <- fda::Data2fd(y = obs_data,argvals = seq(0,1,length.out=length(obs_data)))

  BTC_Monthly[,i] <- as.numeric(fda::eval.fd(tmp,seq(0,1,length.out=30)))
}

res <- detect_changepoint_final_Tn(X = BTC_Monthly)
res <- detect_changepoint_final_Tn(X = BTC_Monthly,h = 5)
res$pval
res1 <- detect_changepoint_Welch(X = BTC_Monthly)
res1 <- detect_changepoint_Welch(X = BTC_Monthly,h = 5)
res1$detected
res2 <- detect_changepoint_bootstrap(X = BTC_Monthly,statistic = 'Tn')
res2 <- detect_changepoint_bootstrap(X = BTC_Monthly,statistic = 'Tn',blockSize = 5)
res2$pval

##########################
MSFT <- read.csv("C:/Users/jerem/Downloads/stock/MSFT.csv")
MSFT$Date <- as.POSIXct(MSFT$Local.time,format='%d.%m.%Y %H:%M:%S')
MSFT$YearMonth <- paste0(year(MSFT$Date),'_',month(MSFT$Date))
MSFT$YearMonthDay <- paste0(year(MSFT$Date),'_',month(MSFT$Date),'_',day(MSFT$Date))
MSFT$Day <- day(MSFT$Date)
MSFT$Hour <- hour(MSFT$Date)

MSFT_Daily <- data.frame(matrix(NA,ncol=length(unique(MSFT$YearMonthDay)),
                                 nrow=24))
colnames(MSFT_Daily) <- unique(MSFT$YearMonthDay)
for(i in 1:ncol(MSFT_Daily)){
  obs_data <- MSFT[MSFT$YearMonthDay==colnames(MSFT_Daily)[i],"Close"]
  if(length(obs_data)==24){
    MSFT_Daily[,i] <- obs_data
  } else if(length(obs_data) < 24){
    MSFT_Daily[,i] <- c(obs_data[1:2],mean(obs_data[2],obs_data[3]),obs_data[3:23])
  } else{
    MSFT_Daily[,i] <- c(obs_data[1],mean(obs_data[2:3]),obs_data[4:25])
  }
}

# Clean
MSFT_Daily <- MSFT_Daily[9:16,]  # Remove non-trading hrs
uniquelength <- sapply(MSFT_Daily,function(x) length(unique(x)))
MSFT_Daily <- subset(MSFT_Daily, select=uniquelength>1) # Remove non-trading days

bs <- binary_segmentation(MSFT_Daily,'Tn','Sim') # all

#

msft_cidr <- MSFT_Daily
for(i in 1:nrow(MSFT_Daily)){
  msft_cidr[i,] <- log(MSFT_Daily[i,]) - log(MSFT_Daily[1, ])
}

bs_cidr <- binary_segmentation(msft_cidr,'Tn','Sim')

#
msft_fd <- fda::Data2fd(1:8, as.matrix(MSFT_Daily))

# Play with more components
nPCs <- 5
msft_fpca <- fda::pca.fd(msft_fd, nharm = nPCs)
msft_fpca_comp <- msft_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(msft_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

msft_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- msft_fpca$harmonics$coefs %*% t(msft_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:8, msft_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')

##########################

nox <- read.csv("C:/Users/jerem/Downloads/nox.csv")[,-1]
nox$date <- as.Date(nox$date,'%Y-%m-%d')
# 2020-2022
nox0 <- nox[year(nox$date) %in% c(2018:2023),]
nox0 <- t(nox0)
colnames(nox0) <- nox0[1,]
nox0 <- as.matrix(nox0[-1,])

use_dat <- apply(nox0,MARGIN = 2, as.numeric)
nox_int <- linear_imputatation(use_dat,use.prev.curve = T)

nox_fd <- fda::Data2fd(1:nrow(nox_int), nox_int)
# Play with more components
nPCs <- 10
nox_fpca <- fda::pca.fd(nox_fd, nharm = nPCs)
nox_fpca_comp <- nox_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(nox_fpca_comp[, i], freq = 7)
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}
nox_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- nox_fpca$harmonics$coefs %*% t(nox_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:nrow(nox_int), nox_fd$basis) %*% orig_coefs

bs_nox_fpca <- binary_segmentation(X = eval_fd_vals,
                                   statistic = 'Tn', method = 'Sim')

plot_fd(use_dat,bs_nox_fpca)

##########################

HPI <- read.csv("C:/Users/jerem/Downloads/TryCE/HPI.csv")
HPI_wide <- HPI %>% pivot_wider(id_cols = Quarter,
                                names_from = Year,
                                values_from = HPI)
HPI_wide <- as.data.frame(HPI_wide[,-c(1,ncol(HPI_wide))])

bs_HPI_fpca <- binary_segmentation(X = HPI_wide,
                                   statistic = 'Tn', method = 'Sim')

##########################

House <- read.csv("C:/Users/jerem/Downloads/TryCE/CSUSHPINSA.csv")
House$DATE <- as.Date(House$DATE,"%m/%d/%Y")
House$Year <- year(House$DATE)
House$Month <- month(House$DATE)
House_wide <- House %>% pivot_wider(id_cols = Month,
                                    names_from = Year,
                                    values_from = CSUSHPINSA)
House_wide <- as.data.frame(House_wide[,-1])

##########################
rates <- read.csv("C:/Users/jerem/Downloads/TryCE/yieldCurves.csv")
rates_prep <- t(rates)
colnames(rates_prep) <- rates_prep[1,]
rates_prep <- rates_prep[-1,rev(1:ncol(rates_prep))]
rates_prep <- apply(rates_prep,MARGIN = 2, as.numeric)
rates_prep_int <- linear_imputatation(rates_prep,use.prev.curve = T)

bs_rates_fpca <- binary_segmentation(X = rates_prep_int,
                                   statistic = 'Tn', method = 'Sim')

rates_prep_int_23 <- rates_prep_int[,
  as.Date(colnames(rates_prep_int),'%m/%d/%Y')>=as.Date("2023-01-01")]

bs_rates23_fpca <- binary_segmentation(
  X = rates_prep_int_23, statistic = 'Tn', method = 'Sim')
plot_fd(rates_prep_int_23,bs_rates23_fpca)



rates_fd <- fda::Data2fd(c(1,2,3,4,6,12,24,36,60,84,120,240,360),
                         as.matrix(rates_prep_int_23))

# Play with more components
nPCs <- 5
rates_fpca <- fda::pca.fd(rates_fd, nharm = nPCs)
rates_fpca_comp <- rates_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(rates_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

rates_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- rates_fpca$harmonics$coefs %*% t(rates_fpca_forecast)
eval_fd_vals <- fda::eval.basis(c(1,2,3,4,6,12,24,36,60,84,120,240,360),
                                rates_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
##########################

prec <- read.csv("C:/Users/jerem/Downloads/TryCE/LAPrep.csv")
prec$DATE <- as.Date(prec$DATE,"%m/%d/%Y")
prec$YearMonth <- paste0(year(prec$DATE),'_',month(prec$DATE))
prec$Day <- day(prec$DATE)
prec_wide <- data.frame(matrix(ncol=length(unique(prec$YearMonth)),
                               nrow=30))
colnames(prec_wide) <- unique(prec$YearMonth)

for(i in 1:ncol(prec_wide)){
  tmp <- prec[prec$YearMonth==colnames(prec_wide)[i],"PRCP"]
  obs_data <- #suppressWarnings(
    fda::Data2fd(
      argvals = seq(0,1,length.out=length(tmp)),
      y=tmp,
      basisobj = fda::create.bspline.basis(c(0,1),20)
    )
  #)

  prec_wide[,i] <- pmax(0,as.numeric(fda::eval.fd(seq(0,1,length.out=30),obs_data)))

}

# prec_wide <- prec %>% pivot_wider(id_cols = Day,
#                                     names_from = YearMonth,
#                                     values_from = PRCP)
# prec_wide <- prec_wide[order(prec_wide$Day),-1]
# prec_wide[is.na(prec_wide)] <- 0


bs_prec_fpca <- binary_segmentation(
  X = prec_wide, statistic = 'Tn', method = 'Sim')

##########################

prec <- read.csv("C:/Users/jerem/Downloads/TryCE/HiPrec.csv")
prec$DATE <- as.Date(prec$DATE,"%Y-%m-%d")
prec$YearMonth <- paste0(year(prec$DATE),'_',month(prec$DATE))
prec$Day <- day(prec$DATE)
prec_wide <- data.frame(matrix(ncol=length(unique(prec$YearMonth)),
                               nrow=30))
colnames(prec_wide) <- unique(prec$YearMonth)

for(i in 1:ncol(prec_wide)){
  tmp <- prec[prec$YearMonth==colnames(prec_wide)[i],c('PRCP','YearMonth','Day')]
  tmp <- na.omit(tmp)

  obs_data <- #suppressWarnings(
    fda::Data2fd(
      argvals = tmp$Day/max(tmp$Day),
      y=tmp$PRCP,
      basisobj = fda::create.bspline.basis(c(0,1),nrow(tmp))
    )
  #)

  prec_wide[,i] <- pmax(0,as.numeric(fda::eval.fd(seq(0,1,length.out=30),obs_data)))

}

# prec_wide <- prec %>% pivot_wider(id_cols = Day,
#                                     names_from = YearMonth,
#                                     values_from = PRCP)
# prec_wide <- prec_wide[order(prec_wide$Day),-1]
# prec_wide[is.na(prec_wide)] <- 0


bs_HiPrec_fpca <- binary_segmentation(
  X = prec_wide, statistic = 'Tn', method = 'Sim')

##########################

prec <- read.csv("C:/Users/jerem/Downloads/TryCE/HiPrec.csv")
prec$DATE <- as.Date(prec$DATE,"%Y-%m-%d")
prec$Year <- year(prec$DATE)
prec$Month <- month(prec$DATE)
prec$Day <- day(prec$DATE)
prec_wide <- data.frame(matrix(ncol=length(unique(prec$Year)),
                               nrow=360))
colnames(prec_wide) <- unique(prec$Year)

for(i in 1:ncol(prec_wide)){
  tmp <- prec[prec$Year==colnames(prec_wide)[i],c('PRCP','Year','Month','Day')]
  tmp <- na.omit(tmp)
  tmp <- tmp[order(tmp$Year,tmp$Month,tmp$Day),]

  obs_data <- #suppressWarnings(
    fda::Data2fd(
      argvals = tmp$Day/max(tmp$Day),
      y=tmp$PRCP,
      basisobj = fda::create.bspline.basis(c(0,1),nrow(tmp))
    )
  #)

  prec_wide[,i] <- pmax(0,
                        as.numeric(
                          fda::eval.fd(
                            seq(0,1,length.out=360),obs_data)))

}

bs_HiPrec_fpca <- binary_segmentation(
  X = prec_wide, statistic = 'Tn', method = 'Sim')

##########################

prec <- read.csv("C:/Users/jerem/Downloads/TryCE/HiPrec.csv")
prec$DATE <- as.Date(prec$DATE,"%Y-%m-%d")
prec$Year <- year(prec$DATE)
prec$Month <- month(prec$DATE)
prec$Day <- day(prec$DATE)
prec_wide <- data.frame(matrix(ncol=length(unique(prec$Year)),
                               nrow=12))
colnames(prec_wide) <- unique(prec$Year)

for(i in 1:ncol(prec_wide)){
  tmp <- prec[prec$Year==colnames(prec_wide)[i],c('PRCP','Year','Month','Day')]
  tmp <- na.omit(tmp)
  tmp <- tmp[order(tmp$Year,tmp$Month,tmp$Day),]

  val <- rep(NA,12)
  for(j in 1:12){
    val[j] <- mean(tmp[tmp$Month==j,"PRCP"])
    if(is.na(val[j])){
      val[j] <- rowMeans(prec_wide,na.rm = T)[j]
    }
  }


  obs_data <- #suppressWarnings(
    fda::Data2fd(
      argvals = seq(0,1,length.out=12),
      y=val,
      basisobj = fda::create.bspline.basis(c(0,1),12)
    )
  #)

  prec_wide[,i] <- pmax(0,
                        as.numeric(
                          fda::eval.fd(
                            seq(0,1,length.out=12),obs_data)))

}

bs_HiPrec_fpca <- binary_segmentation(
  X = prec_wide, statistic = 'Tn', method = 'Sim')

##########################

tsla <- read.csv("C:/Users/jerem/Downloads/TryCE/stock/tsla.csv")
tsla$Date <- as.POSIXct(tsla$Local.time,format='%d.%m.%Y %H:%M:%S')
tsla$YearMonth <- paste0(year(tsla$Date),'_',month(tsla$Date))
tsla$YearMonthDay <- paste0(year(tsla$Date),'_',month(tsla$Date),'_',day(tsla$Date))
tsla$Day <- day(tsla$Date)
tsla$Hour <- hour(tsla$Date)

tsla_Daily <- data.frame(matrix(NA,
                                ncol=length(unique(tsla$YearMonthDay)),
                                nrow=24))
colnames(tsla_Daily) <- unique(tsla$YearMonthDay)
for(i in 1:ncol(tsla_Daily)){
  obs_data <- tsla[tsla$YearMonthDay==colnames(tsla_Daily)[i],"Close"]
  if(length(obs_data)==24){
    tsla_Daily[,i] <- obs_data
  } else if(length(obs_data) < 24){
    tsla_Daily[,i] <- c(obs_data[1:2],mean(obs_data[2],obs_data[3]),obs_data[3:23])
  } else{
    tsla_Daily[,i] <- c(obs_data[1],mean(obs_data[2:3]),obs_data[4:25])
  }
}

# Clean
tsla_Daily <- tsla_Daily[9:16,]  # Remove non-trading hrs
uniquelength <- sapply(tsla_Daily,function(x) length(unique(x)))
tsla_Daily <- subset(tsla_Daily, select=uniquelength>1) # Remove non-trading days

bs <- binary_segmentation(X=tsla_Daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(tsla_Daily,bs)

#

tsla_cidr <- tsla_Daily
for(i in 1:nrow(tsla_Daily)){
  tsla_cidr[i,] <- log(tsla_Daily[i,]) - log(tsla_Daily[1, ])
}

bs_cidr <- binary_segmentation(tsla_cidr,'Tn','Sim')

#
tsla_fd <- fda::Data2fd(1:8, as.matrix(tsla_Daily))

# Play with more components
nPCs <- 5
tsla_fpca <- fda::pca.fd(tsla_fd, nharm = nPCs)
tsla_fpca_comp <- tsla_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(tsla_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

tsla_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- tsla_fpca$harmonics$coefs %*% t(tsla_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:8, tsla_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
plot_fd(tsla_Daily,bs_fpca)

##########################

# AAPL
tmp <- read.csv("C:/Users/jerem/Downloads/TryCE/stock/AAPL.csv")
tmp1 <- read.csv("C:/Users/jerem/Downloads/TryCE/stock/AAPL0.csv")
tmp <- rbind(tmp1,tmp)
tmp$Date <- as.POSIXct(tmp$Local.time,format='%d.%m.%Y %H:%M:%S')
tmp$YearMonth <- paste0(year(tmp$Date),'_',month(tmp$Date))
tmp$YearMonthDay <- paste0(year(tmp$Date),'_',month(tmp$Date),'_',day(tmp$Date))
tmp$Day <- day(tmp$Date)
tmp$Hour <- hour(tmp$Date)
tmp <- tmp[tmp$Date>=as.Date('2021-01-01'),]

tmp_Daily <- data.frame(matrix(NA,
                                ncol=length(unique(tmp$YearMonthDay)),
                                nrow=24))
colnames(tmp_Daily) <- unique(tmp$YearMonthDay)
for(i in 1:ncol(tmp_Daily)){
  obs_data <- tmp[tmp$YearMonthDay==colnames(tmp_Daily)[i],"Close"]
  if(length(obs_data)==24){
    tmp_Daily[,i] <- obs_data
  } else if(length(obs_data) < 24){
    tmp_Daily[,i] <- c(obs_data[1:2],mean(obs_data[2],obs_data[3]),obs_data[3:23])
  } else{
    tmp_Daily[,i] <- c(obs_data[1],mean(obs_data[2:3]),obs_data[4:25])
  }
}

# Clean
tmp_Daily <- tmp_Daily[9:16,]  # Remove non-trading hrs
uniquelength <- sapply(tmp_Daily,function(x) length(unique(x)))
tmp_Daily <- subset(tmp_Daily, select=uniquelength>1) # Remove non-trading days
dates <- colnames(tmp_Daily)
tmp_Daily <- linear_imputatation(tmp_Daily,use.prev.curve = T)
colnames(tmp_Daily) <- dates

plot_fd(tmp_Daily)


# bs <- binary_segmentation(X=tmp_Daily,
#                           statistic = 'Tn',
#                           method = 'Sim') # all
#
# plot_fd(tmp_Daily,bs)

#

tmp_cidr <- tmp_Daily
for(i in 1:nrow(tmp_Daily)){
  tmp_cidr[i,] <- log(tmp_Daily[i,]) - log(tmp_Daily[1, ])
}

bs_cidr <- binary_segmentation(tmp_cidr,'Tn','Sim')

#
tmp_fd <- fda::Data2fd(1:8, as.matrix(tmp_Daily))

# Play with more components
nPCs <- 5
tmp_fpca <- fda::pca.fd(tmp_fd, nharm = nPCs)
tmp_fpca_comp <- tmp_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(tmp_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

tmp_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- tmp_fpca$harmonics$coefs %*% t(tmp_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:8, tmp_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
plot_fd(tmp_Daily,bs_fpca)

##########################
OJ <- read.csv("C:/Users/jerem/Downloads/TryCE/OJ.csv")
OJ$Date <- as.POSIXct(OJ$Local.time,format='%d.%m.%Y %H:%M:%S')
OJ$YearMonth <- paste0(year(OJ$Date),'_',month(OJ$Date))
OJ$YearMonthDay <- paste0(year(OJ$Date),'_',month(OJ$Date),'_',day(OJ$Date))
OJ$Day <- day(OJ$Date)
OJ$Hour <- hour(OJ$Date)

OJ_Daily <- data.frame(matrix(NA,
                               ncol=length(unique(OJ$YearMonthDay)),
                               nrow=24))
colnames(OJ_Daily) <- unique(OJ$YearMonthDay)
for(i in 1:ncol(OJ_Daily)){
  obs_data <- OJ[OJ$YearMonthDay==colnames(OJ_Daily)[i],"Close"]
  if(length(obs_data)==24){
    OJ_Daily[,i] <- obs_data
  } else if(length(obs_data) < 24){
    OJ_Daily[,i] <- c(obs_data[1:2],mean(obs_data[2],obs_data[3]),obs_data[3:23])
  } else{
    OJ_Daily[,i] <- c(obs_data[1],mean(obs_data[2:3]),obs_data[4:25])
  }
}

# Clean
OJ_Daily <- OJ_Daily[9:16,]  # Remove non-trading hrs
uniquelength <- sapply(OJ_Daily,function(x) length(unique(x)))
OJ_Daily <- subset(OJ_Daily, select=uniquelength>1) # Remove non-trading days
dates <- colnames(OJ_Daily)
OJ_Daily <- linear_imputatation(OJ_Daily,use.prev.curve = T)
colnames(OJ_Daily) <- dates

plot_fd(OJ_Daily)

##########################
XAUUSD <- read.csv("C:/Users/jerem/Downloads/TryCE/XAUUSD.csv")
XAUUSD$Date <- as.POSIXct(XAUUSD$Local.time,format='%d.%m.%Y %H:%M:%S')
XAUUSD$YearMonth <- paste0(year(XAUUSD$Date),'_',month(XAUUSD$Date))
XAUUSD$YearMonthDay <- paste0(year(XAUUSD$Date),'_',month(XAUUSD$Date),'_',day(XAUUSD$Date))
XAUUSD$Day <- day(XAUUSD$Date)
XAUUSD$Hour <- hour(XAUUSD$Date)

XAUUSD_Daily <- data.frame(matrix(NA,
                              ncol=length(unique(XAUUSD$YearMonthDay)),
                              nrow=24))
colnames(XAUUSD_Daily) <- unique(XAUUSD$YearMonthDay)
for(i in 1:ncol(XAUUSD_Daily)){
  obs_data <- XAUUSD[XAUUSD$YearMonthDay==colnames(XAUUSD_Daily)[i],"Close"]
  if(length(obs_data)==24){
    XAUUSD_Daily[,i] <- obs_data
  } else if(length(obs_data) < 24){
    XAUUSD_Daily[,i] <- c(obs_data[1:2],mean(obs_data[2],obs_data[3]),obs_data[3:23])
  } else{
    XAUUSD_Daily[,i] <- c(obs_data[1],mean(obs_data[2:3]),obs_data[4:25])
  }
}

# Clean
XAUUSD_Daily <- XAUUSD_Daily[9:16,]  # Remove non-trading hrs
uniquelength <- sapply(XAUUSD_Daily,function(x) length(unique(x)))
XAUUSD_Daily <- subset(XAUUSD_Daily, select=uniquelength>1) # Remove non-trading days
dates <- colnames(XAUUSD_Daily)
XAUUSD_Daily <- linear_imputatation(XAUUSD_Daily,use.prev.curve = T)
colnames(XAUUSD_Daily) <- dates

plot_fd(XAUUSD_Daily)

bs <- binary_segmentation(X=XAUUSD_Daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(XAUUSD_Daily,bs)

# XAUUSD_cidr <- XAUUSD_Daily
# for(i in 1:nrow(XAUUSD_Daily)){
#   XAUUSD_cidr[i,] <- log(XAUUSD_Daily[i,]) - log(XAUUSD_Daily[1, ])
# }
#
# bs_cidr <- binary_segmentation(XAUUSD_cidr,'Tn','Sim')

#
XAUUSD_fd <- fda::Data2fd(1:8, as.matrix(XAUUSD_Daily))

# Play with more components
nPCs <- 5
XAUUSD_fpca <- fda::pca.fd(XAUUSD_fd, nharm = nPCs)
XAUUSD_fpca_comp <- XAUUSD_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(XAUUSD_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

XAUUSD_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- XAUUSD_fpca$harmonics$coefs %*% t(XAUUSD_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:8, XAUUSD_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
plot_fd(XAUUSD_Daily,bs_fpca)

##########################
covid <- read.csv("C:/Users/jerem/Downloads/TryCE/covid.csv")
covid$date <- as.Date(covid$date)
covid$YearMonth <- paste0(year(covid$date),'_',month(covid$date))
covid$YearMonthDay <- paste0(year(covid$date),'_',
                             month(covid$date),'_',day(covid$date))
covid$Day <- day(covid$date)

covid_Daily <- data.frame(matrix(NA,
                                  ncol=length(unique(covid$YearMonth)),
                                  nrow=30))
colnames(covid_Daily) <- unique(covid$YearMonth)
for(i in 1:ncol(covid_Daily)){
  tmp <- covid[covid$YearMonth==colnames(covid_Daily)[i],c('hospitalCases','Day')]
  tmp <- na.omit(tmp)
  tmp <- tmp[order(tmp$Day),]

  if(max(tmp$Day)>=28){
    obs_data <- #suppressWarnings(
      fda::Data2fd(
        argvals = seq(0,1,length.out=max(tmp$Day)),
        y=tmp$hospitalCases,
        basisobj = fda::create.bspline.basis(c(0,1),max(tmp$Day))
      )
    #)

    covid_Daily[,i] <- pmax(0,
                          as.numeric(
                            fda::eval.fd(
                              seq(0,1,length.out=30),obs_data)))
  } else{
    covid_Daily[,i] <- NA
  }
}
covid_Daily <- covid_Daily[,-1]

plot_fd(covid_Daily)

##########################
power <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/",
                         "household_power_consumption.txt"),sep = ';')
power$fullTime <- as.POSIXct(paste(power$Date,power$Time),format="%d/%m/%Y %H:%M:%S")
power$YearMonthDay <- paste0(year(power$fullTime),'_',
                             month(power$fullTime),'_',day(power$fullTime))
power$Min <- minute(power$fullTime)

power_Daily <- data.frame(matrix(NA,
                                  ncol=length(unique(power$YearMonthDay)),
                                  nrow=60*24))
colnames(power_Daily) <- unique(power$YearMonthDay)
for(i in 1:ncol(power_Daily)){
  obs_data <- power[power$YearMonthDay==colnames(power_Daily)[i],
                    c('Global_active_power','fullTime')]
  obs_data <- obs_data[order(obs_data$fullTime),]
  if(nrow(obs_data)!=1440){
    next
  } else{
    power_Daily[,i] <- suppressWarnings(
      as.numeric(obs_data$Global_active_power)
    )
  }
}
power_Daily <- power_Daily[,-1]
power_Daily <- linear_imputatation(power_Daily,use.prev.curve = T)

plot_fd(power_Daily)

bs <- binary_segmentation(X=power_Daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(power_Daily,bs)

#
power_fd <- fda::Data2fd(1:nrow(power_Daily), as.matrix(power_Daily))

# Play with more components
nPCs <- 5
power_fpca <- fda::pca.fd(power_fd, nharm = nPCs)
power_fpca_comp <- power_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(power_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

power_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- power_fpca$harmonics$coefs %*% t(power_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:8, power_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
plot_fd(power_Daily,bs_fpca)

##########################
Jan2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/January_2023.csv", na.strings="na")
Feb2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/February_2023.csv", na.strings="na")
Mar2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/March_2023.csv", na.strings="na")
Apr2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/April_2023.csv", na.strings="na")
May2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/May_2023.csv", na.strings="na")
Jun2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/June_2023.csv", na.strings="na")
Jul2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/July_2023.csv", na.strings="na")
Aug2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/August_2023.csv", na.strings="na")
Sep2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/September_2023.csv", na.strings="na")
Oct2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/October_2023.csv", na.strings="na")
Nov2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/November_2023.csv", na.strings="na")
Dec2023 <- read.csv("C:/Users/jerem/Downloads/TryCE/peds/December_2023.csv", na.strings="na")

peds <- plyr::rbind.fill(Jan2023,Feb2023,Mar2023,Apr2023,May2023,Jun2023,
                  Jul2023,Aug2023,Sep2023,Oct2023,Nov2023,Dec2023)
dropCols <- c()
for(i in 3:ncol(peds)){
  peds[,i] <- suppressWarnings(as.numeric(peds[,i]))
  if(sum(is.na(peds[,i]))>0)
    dropCols <- c(dropCols,i)
}
peds <- peds[,-dropCols]

peds$total <- 0
for(i in 1:nrow(peds)){
  peds$total[i] <- sum(peds[i,-c(1:2,ncol(peds))])
}
peds$Date <- as.Date(peds$Date,'%d/%m/%Y')
peds$MonthDay <- paste0(month(peds$Date),'_',day(peds$Date))
peds$MonthDayHour <- paste0(month(peds$Date),'_',day(peds$Date),'_',peds$Hour)
#peds <- peds[,c(ncol(peds),2,ncol(peds)-1)]
peds_daily <- peds %>% pivot_wider(id_cols = Hour,
                                  names_from = MonthDay,
                                  values_from = Bourke.Street.Mall..South.)
peds_daily <- as.data.frame(peds_daily[,-1])

plot_fd(peds_daily)

bs <- binary_segmentation(X=peds_daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(peds_daily,bs)

##########################

river <- as.numeric( read.csv("C:/Users/jerem/Downloads/TryCE/river.txt",
                  header=FALSE,sep = ',') )
river <- data.frame('river'=river)
river$Date <- as.Date(as.Date('1915-01-01'):as.Date('1979-12-31'))
river$YearMonth <- paste0(year(river$Date),'_',month(river$Date))
river$Day <- day(river$Date)

river_Daily <- data.frame(matrix(nrow=30,ncol=length(unique(river$YearMonth))))
colnames(river_Daily) <- unique(river$YearMonth)

for(i in 1:ncol(river_Daily)){
  tmp <- river[river$YearMonth==colnames(river_Daily)[i],c('river','Day')]
  tmp <- na.omit(tmp)
  tmp <- tmp[order(tmp$Day),]

  obs_data <- #suppressWarnings(
    fda::Data2fd(
      argvals = seq(0,1,length.out=nrow(tmp)),
      y=tmp$river,
      basisobj = fda::create.bspline.basis(c(0,1),nrow(tmp))
    )
  #)

  river_Daily[,i] <- pmax(0,
                          as.numeric(
                            fda::eval.fd(
                              seq(0,1,length.out=30),obs_data)))

}

plot_fd(river_Daily)

bs <- binary_segmentation(X=river_Daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(river_Daily,bs)

#
river_fd <- fda::Data2fd(1:30, as.matrix(river_Daily))

# Play with more components
nPCs <- 5
river_fpca <- fda::pca.fd(river_fd, nharm = nPCs)
river_fpca_comp <- river_fpca$scores

## Forecast Each
ts_dat <- list()
comps <- list()
comps_resids <- list()
for (i in 1:nPCs) {
  ts_dat[[i]] <- ts(river_fpca_comp[, i], freq = 7)
  # comps[[i]] <- forecast::ets(ts_dat[[i]])
  # comps_resids[[i]] <- resid(comps[[i]])
  comps[[i]] <- forecast::auto.arima(ts_dat[[i]])
  comps_resids[[i]] <- resid(comps[[i]])
}

river_fpca_forecast <- do.call(cbind, comps_resids)

# Revert Back to fd
orig_coefs <- river_fpca$harmonics$coefs %*% t(river_fpca_forecast)
eval_fd_vals <- fda::eval.basis(1:30, river_fd$basis) %*% orig_coefs

bs_fpca <- binary_segmentation(eval_fd_vals,'Tn','Sim')
plot_fd(river_Daily,bs_fpca)

##########################
# 10.5281/zenodo.3889785
traffic <- read.csv("C:/Users/jerem/Downloads/TryCE/traffic.txt",
                              header=FALSE,sep = ',')
traffic_time <- expand.grid(as.Date(as.Date('2015-01-01'):as.Date('2016-12-31')),0:23)
colnames(traffic_time) <- c('Date','Hour')
traffic_time <- traffic_time[order(traffic_time$Date,traffic_time$Hour),]
traffic_time <- cbind(traffic_time,traffic)
traffic_time$YearMonthDay <- paste0(year(traffic_time$Date),'_'
                                    ,month(traffic_time$Date),'_',
                                    day(traffic_time$Date))
traffic_time$total <- rowSums(traffic_time[,-c(1,2,ncol(traffic_time))])

traffic_time_daily <- traffic_time %>% pivot_wider(id_cols = Hour,
                                                   names_from = YearMonthDay,
                                                   values_from = total)
traffic_time_daily <- as.data.frame(traffic_time_daily[,-1])
#traffic_time_daily <- linear_imputatation(traffic_time_daily,use.prev.curve = T)
#plot_fd(traffic_time_daily)

bs <- binary_segmentation(X=traffic_time_daily,
                          statistic = 'Tn',
                          method = 'Sim') # all
plot_fd(traffic_time_daily,bs)



##########################
power_data <- readRDS('C:/Users/jerem/Downloads/TryCE/anonymous_public_pv_power_data.rds')
power_data1 <- power_data[power_data$unit==1,]

power_data1$MonDayHr <- paste0(month(power_data1$utc), '_',
                                  day(power_data1$utc),'_',
                                  hour(power_data1$utc))
power_data1$Min <- minute(power_data1$utc)

power_hourly <- data.frame(matrix(ncol=length(unique(power_data1$MonDayHr)),
                                  nrow=60))
colnames(power_hourly) <- unique(power_data1$MonDayHr)

for(i in 1:ncol(power_hourly)){
  tmp <- power_data1[power_data1$MonDayHr==colnames(power_hourly)[i],
                     c('max','Min')]
  tmp1 <- data.frame('Min'=0:59,'max'=NA)

  for(j in 0:59){
    tmp2 <- mean(as.numeric(tmp[tmp$Min==j,][['max']]),na.rm = T)
    tmp1[j,'max'] <- ifelse(is.nan(tmp2),NA,tmp2)
  }

  power_hourly[,i] <- as.numeric(tmp1$max)
}

plot_fd(power_hourly)
power_hourly <- linear_imputatation(power_hourly,use.prev.curve = T)

# power_hourly <- power_data1 %>% pivot_wider(names_from = MonDayHr,
#                                             id_cols = Min,
#                                             values_from = max)
# power_hourly <- as.data.frame(power_hourly[,-1])

bs <- binary_segmentation(X=power_hourly,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(river_Daily,bs)

##########################
power_data <- readRDS('C:/Users/jerem/Downloads/TryCE/anonymous_public_pv_power_data.rds')
power_data1 <- power_data[power_data$unit==1,]

power_data1$MonDay <- paste0(month(power_data1$utc), '_',
                               day(power_data1$utc))
power_data1$Hour <- hour(power_data1$utc)
power_data1$Min <- minute(power_data1$utc)

power_Daily <- data.frame(matrix(ncol=length(unique(power_data1$MonDay)),
                                  nrow=60*24))
colnames(power_Daily) <- unique(power_data1$MonDay)

for(i in 1:ncol(power_Daily)){
  tmp <- power_data1[power_data1$MonDay==colnames(power_Daily)[i],
                     c('max','Hour','Min')]
  tmp1 <- expand.grid(0:23,0:59)
  colnames(tmp1) <- c('Hour','Min')
  tmp1 <- tmp1[order(tmp1$Hour,tmp1$Min),]
  rownames(tmp1) <- 1:nrow(tmp1)
  tmp1$max <- NA

  for(h in 0:23){
    for(m in 0:59){
      tmp2 <- mean(as.numeric(tmp[tmp$Min==m &
                                    tmp$Hour==h,][['max']]),na.rm = T)
      tmp1[tmp1$Hour==h & tmp1$Min==m,'max'] <-
        ifelse(is.nan(tmp2),NA,tmp2)
    }
  }

  power_Daily[,i] <- as.numeric(tmp1$max)
}

plot_fd(power_Daily)
power_Daily <- linear_imputatation(power_Daily,use.prev.curve = T)

# power_hourly <- power_data1 %>% pivot_wider(names_from = MonDayHr,
#                                             id_cols = Min,
#                                             values_from = max)
# power_hourly <- as.data.frame(power_hourly[,-1])

bs <- binary_segmentation(X=power_Daily,
                          statistic = 'Tn',
                          method = 'Sim') # all

plot_fd(power_Daily,bs)

##########################

peds <- data.frame()

for(i in c(2019,2020,2021)){

  Jan <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         'January_',i,".csv"), na.strings="na")
  Feb <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "February_",i,".csv"), na.strings="na")
  Mar <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "March_",i,".csv"), na.strings="na")
  Apr <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "April_",i,".csv"), na.strings="na")
  May <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "May_",i,".csv"), na.strings="na")
  Jun <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "June_",i,".csv"), na.strings="na")
  Jul <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "July_",i,".csv"), na.strings="na")
  Aug <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "August_",i,".csv"), na.strings="na")
  Sep <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "September_",i,".csv"), na.strings="na")
  Oct <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "October_",i,".csv"), na.strings="na")
  Nov <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "November_",i,".csv"), na.strings="na")
  Dec <- read.csv(paste0("C:/Users/jerem/Downloads/TryCE/peds/big/",
                         "December_",i,".csv"), na.strings="na")
  if(i==2020){
    Nov$Date <- paste0(Nov$Date,'20')
    Dec$Date <- paste0(Dec$Date,'20')
  }
  if(i==2021){
    Jul$Date <- stringr::str_replace_all(
      as.character(as.Date(Jul$Date,format='%d/%m/%y')),'-','/')
    Jul$Date <- paste0(substr(Jul$Date,9,10),'/07/2021')
  }

  peds <- plyr::rbind.fill(peds,Jan,Feb,Mar,Apr,May,Jun,
                           Jul,Aug,Sep,Oct,Nov,Dec)
}
dropCols <- c()
for(i in 3:ncol(peds)){
  peds[,i] <- suppressWarnings(as.numeric(peds[,i]))
  peds[,i]<-ifelse(peds[,i]==-1,NA,peds[,i])
  if(sum(is.na(peds[,i]))>500)
    dropCols <- c(dropCols,i)
}
peds <- peds[,-dropCols]

peds$total <- 0
for(i in 1:nrow(peds)){
  peds$total[i] <- sum(peds[i,-c(1:2,ncol(peds))],na.rm = T)
}
peds$Date <- as.Date(peds$Date,'%d/%m/%Y')
peds$YearMonthDay <- paste0(year(peds$Date),'_',month(peds$Date),'_',day(peds$Date))
peds$MonthDayHour <- paste0(month(peds$Date),'_',day(peds$Date),'_',peds$Hour)
#peds <- peds[,c(ncol(peds),2,ncol(peds)-1)]

peds_daily <- peds %>% pivot_wider(id_cols = Hour,
                                   names_from = YearMonthDay,
                                   values_from = total)
peds_daily <- as.data.frame(peds_daily[,-1])
peds_daily <- linear_imputatation(peds_daily,use.prev.curve = T)

#plot_fd(peds_daily)
bs <- binary_segmentation(X=peds_daily,
                          statistic = 'Tn',
                          method = 'Sim') # all
plot_fd(peds_daily,bs)
# Bourke.Street.Mall..South.
# Southern.Cross.Station
# Chinatown.Lt.Bourke.St..South.
# Lonsdale.St..South.
# QVM.Therry.St..South.
#! Elizabeth.St.Lonsdale.St..South.

##########################

