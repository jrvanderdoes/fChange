library(tidyverse)

# faithful_data <- read.delim("C:/Users/jerem/Downloads/Old_Faithful_eruptions.tsv")
faithful_data <- read.delim('C:/Users/jerem/Downloads/geysertimes_eruptions_complete_2024-10-03.tsv')

faithful_data$time <-
  as.Date(as.POSIXct(faithful_data$eruption_time_epoch, origin="1970-01-01"))

tmp <- faithful_data[,c("geyser","duration_seconds","time")]
# tmp <- tmp[tmp$duration_seconds>0,]

tmp_unique <- expand.grid(1:12, 1970:2023 )
colnames(tmp_unique) <- c('month', 'year')
tmp_unique$eruptions <- 0
for(i in 1:nrow(tmp_unique)){
  yr <- tmp_unique$year[i]
  mn <- tmp_unique$month[i]
  tmp_unique[i,"eruptions"] <- nrow(tmp[lubridate::year(tmp$time) == yr &
                                          lubridate::month(tmp$time) == mn, ])
}

wide_faithful <- tmp_unique %>%
  pivot_wider(id_cols = month, names_from = year,values_from = eruptions)

plot_fd(funts(wide_faithful), CPs=c(32))
