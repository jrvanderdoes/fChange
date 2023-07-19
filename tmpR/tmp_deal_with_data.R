### Electricity

electricity <- t(fdaACF::elec_prices)

### AustTemp
aust_tmp <- readRDS('C:/Users/Jeremy/Downloads/Australian_Temp.RDS')
aust_tmp <- aust_tmp$Sydney

aust_tmp  <- aust_tmp[,-c(1,2,7,8)]
colnames(aust_tmp) <- c('Year','Month','Day','MinTemp')

aust_tmp1 <- tidyr::pivot_wider(aust_tmp,
                                names_from = c('Year'),
                                values_from = c('MinTemp'))
aust_tmp1 <- aust_tmp1[-366,]
#plot_fd(as.data.frame(aust_tmp1[,-c(1:2)]), 1:365)

aust_tmp <-
  tidyr::`%>%`(
    tidyr::`%>%`(
      tidyr::`%>%`(
        aust_tmp,
        tidyr::unite("date", Year, Month, Day, sep = "-", remove = FALSE)
      ),
      dplyr::mutate(date = lubridate::ymd(date))
    ),
    #dplyr::mutate(doy = lubridate::yday(date))
    dplyr::mutate(doy = (lubridate::yday(date)-1) /
                    (lubridate::yday(lubridate::ymd(paste0(substr(date,1,4),'-12-31')))-1 ) )
  )
aust_tmp  <- aust_tmp[,-c(1,3,4)]
colnames(aust_tmp) <- c('Year','MinTemp','DOY')

aust_tmp <- tidyr::pivot_wider(aust_tmp,
                               names_from = c('Year'),
                               values_from = c('MinTemp'))
aust_tmp <- as.data.frame(aust_tmp)
evalPts <- seq(0,1,0.001)
aust_tmp_evaled <- matrix(nrow=length(evalPts), ncol=ncol(aust_tmp)-1)
for(i in 2:ncol(aust_tmp)){
  tmp_data <- stats::na.omit(aust_tmp[,c(1,i)])
  tmp <- fda::Data2fd(tmp_data[,1], tmp_data[,2],
                      basisobj = fda::create.bspline.basis(nbasis = 21))
  aust_tmp_evaled[,i-1] <- fda::eval.fd(evalPts,tmp)
}
# Drop last column since not complete
aust_tmp1 <- as.data.frame(aust_tmp1[,-c(1:2,ncol(aust_tmp1))])
aust_tmp_evaled <- aust_tmp_evaled[,-ncol(aust_tmp_evaled)]


aust_tmp2 <- linear_imputatation(aust_tmp1)


