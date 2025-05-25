load(file = "./data-raw/electricity.RData")
# colnames(electricity) <-
#   as.character(seq(from = as.Date('01-01-2014','%m-%d-%Y'),
#       to = as.Date('12-31-2014','%m-%d-%Y'),
#       by = "days"))

electricity <- dfts(electricity,name='Spanish Electricity Spot Prices 2014',
                    labels = as.character(seq(from = as.Date('01-01-2014','%m-%d-%Y'),
                                              to = as.Date('12-31-2014','%m-%d-%Y'),
                                              by = "days")),
                    season = 7)

use_data(electricity,overwrite = TRUE)
