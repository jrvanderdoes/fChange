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
                                   values_from = Elizabeth.St.Lonsdale.St..South.)
peds_daily <- as.data.frame(peds_daily[,-1])
peds_daily <- linear_imputatation(peds_daily,use.prev.curve = T)

#plot_fd(peds_daily)
bs <- binary_segmentation(X=peds_daily,
                          statistic = 'Tn',
                          method = 'Sim') # all
plot_fd(peds_daily,bs)

