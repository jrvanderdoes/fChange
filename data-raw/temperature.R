library(tidyverse)
load(file = "./data-raw/Australian_Temp.RData")

temperature <- list()

cities <- names(Australian_Temp)
for (i in 1:length(cities)) {
  data_city <- Australian_Temp[[cities[i]]] %>%
    as.data.frame()
  colnames(data_city) <- c(
    "Product", "Station",
    "Year", "Month", "Day",
    "MinTemp", "DaysAccum",
    "Quality"
  )

  # stations <- unique(data_city$Station)
  # for(j in 1:length(stations)){
  # Drop all but key vars
  data_station <- data_city %>%
    # filter(Station==stations[j]) %>%
    select(-c(
      "Product", "Station",
      "DaysAccum", "Quality"
    ))
  # Organize into fChange shape
  data_station <- data_station %>%
    arrange(Year, Month, Day) %>%
    pivot_wider(names_from = Year, values_from = MinTemp) %>%
    filter(!(Day == 29 & Month == 2)) %>%
    select(-c("Month", "Day")) %>%
    as.data.frame()

  # Fill in missing values
  data_station <- linear_imputatation(
    data = data_station,
    use.prev.curve = TRUE
  )
  # temperature[[paste0(cities[i],"_",stations[j])]] <-
  temperature[[cities[i]]] <-
    data_station
  # }
}

use_data(temperature)
