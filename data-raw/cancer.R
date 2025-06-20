# USA Breast Cancer Data from WHO
load(file = "./data-raw/cancer.RData")

cancer <- cancer %>%
  tidyr::pivot_wider(
    id_cols = `Age Group`, names_from = Year,
    values_from = `Percentage of cause-specific deaths out of total deaths`
  )
cnames <- colnames(cancer)

cancer <- cancer %>% data.frame()
rownames(cancer) <- cancer[, 1]
cancer <- cancer[, -1]
colnames(cancer) <- cnames[-1]

use_data(cancer)
