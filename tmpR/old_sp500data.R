SPYUS20192020 <-
  read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                  "SPYUS20192020.csv"))
SPYUS2021 <-
  read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                  "SPYUS2021.csv"))
SPYUS2022 <-
  read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                  "SPYUS2022.csv"))
SPYUS2023 <-
  read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                  "SPYUS2023.csv"))

SPYUS <- rbind(SPYUS20192020,
               SPYUS2021,
               SPYUS2022,
               SPYUS2023)
SPYUS <- SPYUS[SPYUS$Volume!=0,]
SPYUS$datetime <-
  as.POSIXct(SPYUS$Local.time,
             format="%d.%m.%Y %H:%M:%S")
SPYUS$date <- as.Date(SPYUS$datetime)
SPYUS$hourMin <- paste(hour(SPYUS$datetime),'-',minute(SPYUS$datetime))
SPYUS$hr <- hour(SPYUS$datetime)
SPYUS$min <- minute(SPYUS$datetime)
SPYUS_wide <- SPYUS %>%
  pivot_wider(names_from = date,
              values_from = Close,
              id_cols = c(hr,min) ) %>%
  arrange(.,hr,min)

for(i in 3:ncol(SPYUS_wide)){
  # If any missing
  if(sum(is.na(SPYUS_wide[,i]))>0){
    missing_idx <- which(is.na(SPYUS_wide[,i]))

    nonmissing <- na.omit(SPYUS_wide[,i])[[1]]
    first_val <- nonmissing[1]
    #lastval_val <- nonmissing[length(nonmissing)]

    for(j in missing_idx){
      if(j==1) SPYUS_wide[j,i] <- first_val
      else SPYUS_wide[j,i] <- SPYUS_wide[j-1,i]
    }
  }
}

SPYUS_wide <- as.data.frame(SPYUS_wide[,-c(1:2)])

SPYUS_wide_15min <- SPYUS_wide %>%
  mutate(grp = 1+ (row_number()-1) %/% 15) %>%
  group_by(grp) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  select(-grp) %>%
  as.data.frame()

#############
# compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}
compute_ocidr <- function(dat){
  dat_cidr <- dat
  dat_cidr[,1]  <- 100*(log(dat[,1]) - log(dat[1,1]))
  # Obs
  for(j in 2:ncol(dat_cidr)){
    dat_cidr[,j] <- 100 *
      ( log(dat[,j]) - log(dat[nrow(dat_cidr),j-1]) )

    # for(i in 2:nrow(dat_cidr)){
    #   dat_cidr[i,j] <- 100*
    #     ( log(dat[i,j]) - log(dat[1,j]) )
    # }
  }
  dat_cidr
}

# SPYUS_cidr15 <- compute_cidr(SPYUS_wide_15min)
SPYUS_ocidr15 <- compute_ocidr(SPYUS_wide_15min)
SPYUS_ocidr <- compute_ocidr(SPYUS_wide)

set.seed(25641)
spy_min_changes <- binary_segmentation(SPYUS_ocidr,'Tn','Sim')
plot_save <- plotData(SPYUS_ocidr)
plotly::save_image(plot_save,
                   "C:/Users/jerem/Downloads/sp500_ocidr_min_1923.png",
                   width=1600, height=800)
tmp <- paperPlot(SPYUS_ocidr,spy_min_changes)
plotly::save_image(tmp,
                   "C:/Users/jerem/Downloads/sp500_ocidr_min_changes_1923.png",
                   width=1600, height=800)


# plot_save <- plotData(SPYUS_cidr15)
# plotly::save_image(plot_save,
#                    "C:/Users/jerem/Downloads/sp500_cidr.png",
#                    width=1600, height=800)
plot_save <- plotData(SPYUS_ocidr15)
plotly::save_image(plot_save,
                   "C:/Users/jerem/Downloads/sp500_ocidr.png",
                   width=1600, height=800)


set.seed(42513)
spy_changes <- binary_segmentation(SPYUS_ocidr15,'Tn','Sim')

tmp <- paperPlot(SPYUS_ocidr15,spy_changes)
plotly::save_image(tmp,
                   "C:/Users/jerem/Downloads/sp500_ocidr_changes.png",
                   width=1600, height=800)

set.seed(24132)
spy_changes1 <- binary_segmentation(SPYUS_cidr15,'Tn','Sim')
