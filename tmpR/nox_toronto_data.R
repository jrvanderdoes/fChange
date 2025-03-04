# https://www.epa.gov/outdoor-air-quality-data/download-daily-data
library(tidyverse)

filenames <- list.files("C:\\Users\\jerem\\Downloads\\data",
                        pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
df <- do.call("rbind", ldf)
df_plot <- df[,c("Date","Site.ID","Daily.Max.8.hour.CO.Concentration","DAILY_AQI_VALUE")]
df_plot$Date <- as.Date(df_plot$Date,'%m/%d/%Y')
df_plot$Year <- year(df_plot$Date)
df_plot$Day <- day(df_plot$Date)
df_plot$Month <- month(df_plot$Date)

df1 <- df_plot[df_plot$Site.ID=='490350003',
               c("Day","Month","Year",
                 "Daily.Max.8.hour.CO.Concentration")]
df1 <- df1[!(df1$Day==29 & df1$Month==2),]
df1 <- df1[df1$Year!='1981',]
tmp <- df1 %>% pivot_wider(id_cols = c(Day,Month),
                           names_from = Year,
                           values_from = Daily.Max.8.hour.CO.Concentration) %>%
  as.data.frame
plot_fd(tmp[,-c(1:2)],
        FDReps = as.Date( paste0('01-01-',
                                 colnames(tmp[,-c(1,2)])),
                          '%d-%m-%Y'))

df1 <- df_plot[df_plot$Site.ID=='490350003',
               c("Day","Month","Year",
                 "DAILY_AQI_VALUE")]
df1 <- df1[!(df1$Day==29 & df1$Month==2),]
df1 <- df1[df1$Year!='1981',]
tmp <- df1 %>% pivot_wider(id_cols = c(Day,Month),
                           names_from = Year,
                           values_from = DAILY_AQI_VALUE) %>%
  as.data.frame
plot_fd(tmp[,-c(1:2)],
        FDReps = as.Date( paste0('01-01-',
                                 colnames(tmp[,-c(1,2)])),
                          '%d-%m-%Y'))

tmp_int <- linear_imputatation(tmp[,-c(1:2)])
mn1 <- compute_Mn(tmp_int)
plot_fd(tmp_int,changes = mn1$location)
