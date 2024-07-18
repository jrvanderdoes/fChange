########################
##    DATA SETS
########################

########
## Earthquake Detection

# install.packages('eseis')
library(eseis)
eseis::aux_getFDSNdata(
  duration=60*60*24*30
)

########
## PACKAGE OF FDATA

https://cran.r-project.org/web/packages/fds/fds.pdf

########
## Mortality

https://www.mortality.org/Country/Country?cntr=FRATNP
https://cran.r-project.org/web/packages/demography/demography.pdf


########
## AU Temp

http://www.bom.gov.au/


########
## Backtest data

https://www.backtestmarket.com/


################
#   Setup
api_key <- 'SSlObrYTGC7yC7gfFkO9B1SvnyG9I7O0GzbeUrjI'

library("eia")
#library(EIAapi)
library(tidyverse)

########
#   Propane

data <- data.frame()
tmp <- data.frame('fake'=rep(NA,5000))
itr <- 0
while(nrow(tmp)==5000){
  itr <- itr + 1
  cat(itr,'(',nrow(data),'), ')
  tmp <- eia_data(
    dir='petroleum/pri/spt/data/',
    data = c('value'),
    facets = NULL,
    freq = 'daily',
    start = NULL,
    end = NULL,
    sort = NULL,
    length = 5000,
    offset = (itr-1)*5000,
    tidy = TRUE,
    check_metadata = FALSE,
    cache = TRUE,
    key = api_key
  )
  data <- rbind(data,tmp)
}
data_prop <- data[data$`product-name`=='Propane',c("period","value")]
data_prop$period <- as.Date(data_prop$period)

data_prop_wide <- data_prop %>%
  arrange(period) %>%
  mutate('yearmonth'=paste0(year(period),'-',month(period)),
         'day'=day(period))
data_prop_wide <- unique(data_prop_wide)
data_prop_wide <- data_prop_wide %>%
  pivot_wider(id_cols = day,
              names_from = yearmonth,
              values_from = value) %>%
  arrange(day) %>%
  as.data.frame()
saveRDS(data_prop_wide,
        'C:\\Users\\jerem\\Downloads\\fChange-NewData\\Gas\\USPropane.rds')

########
#   Electricity - Operations (Net Generation in CA)
data <- data.frame()
tmp <- data.frame('fake'=rep(NA,5000))
itr <- 0
facets <- list(respondent = "CAL",
               type='NG')
while(nrow(tmp)==5000){
  itr <- itr + 1
  cat(itr,'(',nrow(data),'), ')
  tmp <- eia_data(
    dir='electricity/rto/region-data/data/',
    data = c('value'),
    facets = facets,
    freq = NULL,
    start = NULL,
    end = NULL,
    sort = NULL,
    length = 5000,
    offset =(itr-1)*5000,
    tidy = TRUE,
    check_metadata = FALSE,
    cache = TRUE,
    key = api_key
  )
  data <- rbind(data,tmp)
}
data_op <- data[,c("period","value")]
# as.POSIXct(tmp$period, format="%Y-%m-%dT%H")
# data_op$period <- as.POSIXct(data_op$period, format="%Y-%m-%dT%H")

data_op_wide <- data_op %>%
  arrange(period) %>%
  mutate('yearmonthday'=paste0(year(period),'-',month(period),'-',day(period)),
         'hour'=as.numeric(substr(period,12,14)))
data_op_wide <- unique(data_op_wide)
data_op_wide <- data_op_wide %>%
  pivot_wider(id_cols = hour,
              names_from = yearmonthday,
              values_from = value) %>%
  arrange(hour) %>%
  as.data.frame()
saveRDS(data_op_wide,
        'C:\\Users\\jerem\\Downloads\\fChange-NewData\\Electricity\\USA\\NG_cal.rds')

########
#   Electricity - Operations (Fuel Type in CA)
data <- data.frame()
fuelTypes <- c('COL','NG','NUC','OIL','OTH','SUN','WAT','WND')

for(type_idx in 1:length(fuelTypes)){
  itr <- itr + 1
  cat(fuelTypes[type_idx],': ')

  tmp <- data.frame('fake'=rep(NA,5000))
  itr <- 0
  facets <- list(respondent = "CAL",
                 fueltype = fuelTypes[type_idx])
  while(nrow(tmp)==5000){
    itr <- itr + 1
    cat(itr,'(',nrow(data),'), ')
    tmp <- eia_data(
      dir='electricity/rto/fuel-type-data/data/',
      data = c('value'),
      facets = facets,
      freq = NULL,
      start = NULL,
      end = NULL,
      sort = NULL,
      length = 5000,
      offset =(itr-1)*5000,
      tidy = TRUE,
      check_metadata = FALSE,
      cache = TRUE,
      key = api_key
    )
    data <- rbind(data,tmp)
  }
  cat('\n')
}
data_ft <- data[,c("period","type-name","value")]
# as.POSIXct(tmp$period, format="%Y-%m-%dT%H")
# data_op$period <- as.POSIXct(data_op$period, format="%Y-%m-%dT%H")

# data_ft_wide <- data_ft %>%
#   arrange(period) %>%
#   mutate('yearmonthday'=paste0(year(period),'-',month(period),'-',day(period)),
#          'hour'=as.numeric(substr(period,12,14)))
# data_ft_wide <- unique(data_ft_wide)
# data_ft_wide <- data_ft_wide %>%
#   pivot_wider(id_cols = hour,
#               names_from = yearmonthday,
#               values_from = value) %>%
#   arrange(hour) %>%
#   as.data.frame()
saveRDS(data_ft,
        'C:\\Users\\jerem\\Downloads\\fChange-NewData\\Electricity\\USA\\FT_cal.rds')
########
#   Electricity - Operations (Demand in SoCal)
data <- data.frame()
tmp <- data.frame('fake'=rep(NA,5000))
itr <- 0
facets <- list(parent = "CISO",
               subba = "SCE")
while(nrow(tmp)==5000){
  itr <- itr + 1
  cat(itr,'(',nrow(data),'), ')
  tmp <- eia_data(
    dir='electricity/rto/region-sub-ba-data/data/',
    data = c('value'),
    facets = facets,
    freq = NULL,
    start = NULL,
    end = NULL,
    sort = NULL,
    length = 5000,
    offset =(itr-1)*5000,
    tidy = TRUE,
    check_metadata = FALSE,
    cache = TRUE,
    key = api_key
  )
  data <- rbind(data,tmp)
}
data_d <- data[,c("period","value")]
# as.POSIXct(tmp$period, format="%Y-%m-%dT%H")
# data_op$period <- as.POSIXct(data_op$period, format="%Y-%m-%dT%H")

data_d_wide <- data_d %>%
  arrange(period) %>%
  mutate('yearmonthday'=paste0(year(period),'-',month(period),'-',day(period)),
         'hour'=as.numeric(substr(period,12,14)))
data_d_wide <- unique(data_d_wide)
data_d_wide <- data_d_wide %>%
  pivot_wider(id_cols = hour,
              names_from = yearmonthday,
              values_from = value) %>%
  arrange(hour) %>%
  as.data.frame()
saveRDS(data_d_wide,
        'C:\\Users\\jerem\\Downloads\\fChange-NewData\\Electricity\\USA\\Demand_cal.rds')

########
## XXXXXXX
