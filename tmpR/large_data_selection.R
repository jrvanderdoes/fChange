
# library(devtools)
# load_all()
#
# library(readxl)
# # https://www.eia.gov/petroleum/gasdiesel
# USgas <- read_excel("C:/Users/jerem/Downloads/USgas.xlsx")
# USgas_wide <- USgas %>%
#   pivot_wider(id_cols = QuarterWeek, names_from = YearQuarter,
#               values_from = `Weekly U.S. Regular Reformulated Retail Gasoline Prices  (Dollars per Gallon)`)
# X <- funts(USgas_wide[-14,-1])
# X <- impute(X,method='linear')
# plot_fd(diff(X))
#
# acf(diff(X))
#
# res2 <- eigen_change(diff(X)$data,d = 2,test = 'ind')
#
# plot_fd(diff(X),CPs = res2$change)


###################
library(devtools)
load_all()

library(readxl)
WHO_bcanc <- read_excel("C:/Users/jerem/Downloads/WHO_bcanc.xlsx",
                       sheet = "USA")
# WHO_bcanc <- WHO_bcanc[WHO_bcanc$`Age Group`!+'[85+]',]
WHO_bcanc_wide <- WHO_bcanc %>%
  pivot_wider(id_cols = `Age Group`, names_from = Year,
              values_from = `Percentage of cause-specific deaths out of total deaths`)
X <- funts(WHO_bcanc_wide[,-1])
# summary(X)
plot_fd(X)

res0 <- binary_segmentation(X,statistic = 'Tn',method = 'Approx')
plot_fd(X,CPs=res0[[1]])

plot_fd(diff(X))
res1 <- binary_segmentation(diff(X),statistic = 'Tn',method = 'Approx')
plot_fd(diff(X),CPs=res1[[1]])

# tmp = pcaExploration(X, order=10)
# res3 <- binary_segmentation(tmp$residuals,statistic = 'Tn',method = 'Approx')
# plot_fd(X,CPs=res3[[1]])
# plot_fd(tmp$residuals,CPs=res3[[1]])


# ###################
#
# library(devtools)
# load_all()
#
# library(readxl)
# WHO_card <- read_excel("C:/Users/jerem/Downloads/WHO_card.xlsx",
#                        sheet = "USA_all")
# WHO_card_wide <- WHO_card %>%
#   pivot_wider(id_cols = `Age Group`, names_from = Year, values_from = Number)
# X <- funts(WHO_card_wide[,-1],
#            intraobs = c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85))
# # summary(X)
# plot_fd(X)
#
# res <- mean_change(X)
# res$pvalue
# plot_fd(diff(X))
# res1 <- characteristic_change_sim(diff(X)$data)
#
# res2 <- binary_segmentation(X,statistic = 'Tn',method = 'Sim')
#
# res3 <- covariance_kernel_change(X$data)
# plot_fd(X,CPs=res3$location)
#
# ###################
#
# library(devtools)
# load_all()
#
# library(readxl)
# cancer_au <- read_excel("C:/Users/jerem/Downloads/cancer_au.xlsx",
#                         sheet = "Use")
# cancer_wide <- cancer_au %>%
#   pivot_wider(id_cols = `Age group (years)`, names_from = Year, values_from = Count)
# X <- funts(cancer_wide[,-1])
# #summary(X)
# plot_fd(X)
#
# res <- mean_change(X)
# res$pvalue
# plot_fd(X,CPs = 12)
# plot_fd(diff(X))
# res1 <- characteristic_change_sim(diff(X)$data)
#
# res2 <- binary_segmentation(X,statistic = 'Tn',method = 'Sim')
#
# res3 <- covariance_kernel_change(X$data)
# plot_fd(X,CPs=res3$location)

###################

library(devtools)
load_all()
X <- funts(electricity)
summary(X)


tmp = pcaExploration(X,order=3)
res3 <- binary_segmentation(tmp$residuals,statistic = 'Tn',method = 'Approx')
plot_fd(X,CPs=res3[[1]])
plot_fd(tmp$residuals,CPs=res3[[1]])

###################
# https://home.treasury.gov/interest-rates-data-csv-archive
UScurveyields <- read.csv("C:/Users/jerem/Downloads/UScurveyields.csv")
colnames(UScurveyields) <- c('Date','M1','M2','M3','M4','M6',
                             'M12','M24','M36','M60','M84',
                             'M120','M240','M360')
UScurveyields$Date <- as.Date(UScurveyields$Date,'%m/%d/%y')

library(devtools)
load_all()
X <- funts(t(UScurveyields[,-1]),labels = UScurveyields[,1],
           intraobs = c(1,2,3,4,6,12,24,36,60,84,120,240,360))
X <- funts(X$data[,colSums(is.na(X$data))<nrow(X$data)],
                    labels = X$labels[colSums(is.na(X$data))>0],
                    intraobs = X$intraobs)
# mean_change(X)
X <- impute(X,method='linear')
# res <- binary_segmentation(X = X,'Tn','Sim')
summary(X)


# ###################
# library(tidyverse)
# sweden <- read.csv("C:/Users/jerem/Downloads/SWE/Deaths_1x1.txt", sep="")
# sweden_wide <- sweden %>%
#   pivot_wider(names_from = Year,id_cols = Age,values_from = Total) %>%
#   as.data.frame()
# rownames(sweden_wide) <- sweden_wide[,1]
# sweden_wide <- sweden_wide[,-1]
# swed <- funts(sweden_wide,intraobs = 0:110)
#
#
# result_ce <- binary_segmentation(X = swed,statistic = 'Tn',
#                                  method = 'Approx')
# result <- recursive_segmentation(swed$data)
# plot_fd(swed, result)
# summary(swed,CPs = result)
#
# tmp <- diff(funts(sweden_wide))
# result_ce <- binary_segmentation(X = tmp,statistic = 'Tn', method = 'Sim')
# result <- recursive_segmentation(tmp)
# plot_fd(swed, result)
# summary(swed,CPs = result)

# ###################
# albert_river <- read.delim("C:/Users/jerem/Downloads/albert_river.txt", header=FALSE)
# summary(funts(albert_river))
# plot_fd(t(albert_river),interactive = F)
#
# result_ce <- binary_segmentation(X = t(albert_river),
#                                  statistic = 'Tn',
#                                  method = 'Sim')

# ###################
# mary_river <- read.delim("C:/Users/jerem/Downloads/mary_river.txt", header=FALSE)
# plot_fd(t(mary_river))
#
# result_ce <- binary_segmentation(X = t(mary_river),
#                                  statistic = 'Tn',
#                                  method = 'Sim')
# plot_fd(t(mary_river),CPs=result_ce[[1]])
# ###################
# lockyer_valley <- read.delim("C:/Users/jerem/Downloads/lockyer_valley.txt", header=FALSE)
# plot_fd(t(lockyer_valley))
#
# result_ce <- binary_segmentation(X = t(lockyer_valley),
#                                  statistic = 'Tn',
#                                  method = 'Sim')

# ###################
# temp <- read.delim("C:/Users/jerem/Downloads/temp.txt",
#                    sep='',header=TRUE)
# col_val <- c(temp[,1])
# temp <- t(temp[,-1])
# colnames(temp) <- col_val
# plot_fd(temp)

# ###################
# BL_nom <- read.csv("C:/Users/jerem/Downloads/BL_nom.csv")
# obs <- BL_nom[,1]
# BL_nom <- t(BL_nom[,-1])
# BLC <- funts(BL_nom,intraobs = 1:60/60,labels = as.Date(obs,'%d %b %y'))
# # plot_fd(BLC,interactive = F)
# BLC_impute <- funts(BLC$data[,colSums(is.na(BLC$data))<nrow(BLC$data)],
#                     labels = BLC$labels[colSums(is.na(BLC$data))>0],
#                     intraobs = BLC$intraobs)
# BLC_impute <- impute(BLC_impute, method='mean_obs')
# # plot_fd(BLC_impute,interactive = F)
#
# plot_fd(diff(BLC_impute),interactive = F)
# result_ce <- binary_segmentation(X = diff(BLC_impute),
#                                  statistic = 'Tn',
#                                  method = 'Sim')
# val <- recursive_segmentation(diff(BLC_impute))


# ###################
# ## https://berkeleyearth.org/data/
# BerkleyGlobal <- read.csv("C:/Users/jerem/Downloads/BerkleyGlobal.txt", sep="")
# BerkleyGlobal_wide <- BerkleyGlobal %>%
#   mutate(MonYear=paste0(Month,'-',Year)) %>%
#   pivot_wider(id_cols = Month, names_from = Year, values_from = MonAnomaly)
# X <- funts(BerkleyGlobal_wide[,-1],labels = 1850:2024)
# X <- impute(X,method='linear')
# result_ce <- binary_segmentation(X = X,
#                                  statistic = 'Tn',
#                                  method = 'Sim')
# pt <- plot_fd(X, CPs=result_ce[[1]])
# plotly::save_image(pt,'C:/Users/jerem/Downloads/berkley_global.png')
# val <- recursive_segmentation(X)

# ###################
# ## https://fred.stlouisfed.org/series/DFF
# DFF <- read.csv("C:/Users/jerem/Downloads/DFF.csv")
# DFF <- DFF %>%
#   mutate(DATE=as.Date(DATE)) %>%
#   mutate(day=day(DATE),
#          yearmonth=paste0(year(DATE),'-',month(DATE))) %>%
#   pivot_wider(id_cols = day, names_from = yearmonth,values_from = DFF)
# X <- funts(DFF[,-1],labels = colnames(DFF))
# X <- impute(X,method = 'linear')
# plot_fd(diff(X))
# result_ce1 <- binary_segmentation(X = diff(X),
#                                  statistic = 'Tn',
#                                  method = 'Sim')
# pt1 <- plot_fd(diff(X), CPs=result_ce1[[1]])
# plotly::save_image(pt1,'C:/Users/jerem/Downloads/diff_fedfunds.png')
# pt1.1 <- plot_fd(X, CPs=result_ce1[[1]])
# plotly::save_image(pt1.1,'C:/Users/jerem/Downloads/diff_fedfunds1.png')
