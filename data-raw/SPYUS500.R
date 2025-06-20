# library(tidyverse)
load(file = "./data-raw/SPYUS500.RData")

SPYUS500 <- dfts(SPYUS500,
  name = "S&P500 Index (SPY)",
  labels = as.Date(colnames(SPYUS500))
)

use_data(SPYUS500)
