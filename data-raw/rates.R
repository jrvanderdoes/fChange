# USA Breast Cancer Data from WHO
load(file = "./data-raw/rates.RData")

rates <- dfts(t(rates[nrow(rates):1, -1]),
  labels = rates[nrow(rates):1, 1],
  fparam = c(1, 2, 3, 4, 6, 12, 24, 36, 60, 84, 120, 240, 360),
  name = "US Yield Curves"
)
rates <- dfts(rates$data[, colSums(is.na(rates$data)) < nrow(rates$data)],
  labels = rates$labels[colSums(is.na(rates$data)) < nrow(rates$data)],
  fparam = rates$fparam,
  name = "US Yield Curves"
)

use_data(rates)
