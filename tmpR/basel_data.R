library(tidyverse)

data_path <- 'C:/Users/jerem/OneDrive/Documents/School/Waterloo/Research/GraphCP/Data/Examples/'

tmp <- read.csv2(paste0(data_path,"Basel_Clean.csv"),sep=',')
tmp[["date"]] <- as.Date(tmp$timestamp,"%Y%m%dT%H%M")
tmp[["year"]] <- substr(tmp$timestamp,1,4)
tmp[["day"]] <- yday(tmp$date)
tmp[["time"]] <- substr(tmp$timestamp,10,11)

RH_fd <-
  tmp %>%
  filter(year==2021) %>%
  select('day','time','Relative.Humidity') %>%
  mutate(Relative.Humidity=as.numeric(Relative.Humidity)) %>%
  pivot_wider(id_cols = 'time',
              names_from = 'day',
              values_from = 'Relative.Humidity') %>%
  # filter(!dplyr::row_number() %in% c(391)) %>%
  # mutate(across(`20190102`:`20201231`, as.numeric)) %>%
  select(-c(time)) %>%
  as.data.frame() %>%
  linear_imputatation()


# path2folder <- paste0("C:/Users/jerem/OneDrive/Documents/",
#                       "School/Waterloo/Research/GraphCP/")
# source(file = paste0(path2folder,"covchange.R"))
# source(file = paste0(path2folder,"generalCode.R"))
# source(file = paste0(path2folder,"meanMethod.R"))
# source(file = paste0(path2folder,"tmp_gen_fd.R"))
# source(file = paste0(path2folder,"tmp_plot_fd.R"))

RH_Mn <- compute_Mn(RH_fd)
plot_fd(RH_fd,changes = RH_Mn$location)
plot_fd(RH_fd,val_axis_title = 'Days in 2021')
