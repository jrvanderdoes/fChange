#' Spanish Spot Electricity Data
#'
#' The hourly electricity spot prices from Spain in 2014.
#'
#' @format ## `electricity`
#' A dfts object
#'
#' @source <www.omie.es>
"electricity"

#' Australian Temperature Data
#'
#' A list of daily minimum temperature data from reporting stations in
#'  Australia. Any missing data was previously, linearly, interpolated and leap
#'  days were removed for consistency. Moreover, the data is classified by city
#'  as only one station was observed in each city.
#'
#' @format ## `temperature`
#' A list with 8 elements, one for each city. Each element is a data.frame with
#' 365 rows (one per day), and a variable number of columns related to years of
#' collected data.
#'
#' @source <www.bom.gov.au>
"temperature"

#' S&P 500 Index Data
#'
#' Intraday prices for the S&P500 index (SPY) for 2019 to 2023. Minutely
#'  resolution and daily observations.
#'
#' @format ## `SPY500`
#' A dfts object
"SPYUS500"


#' Breast Cancer
#'
#' Percentage of cause-specific deaths out of total deaths for female breast
#'  cancer in the United States from 1950 to 2021.
#'
#' @format ## `cancer`
#' A data.frame with columns being the year and rows the age groups (5 years)
"cancer"


#' US Yield Curves
#'
#' Yield curves in the US for 1 - 360 month maturity from 1990 to 2023
#'  (removing days without information).
#'
#' @format ## `rates`
#' A dfts object
"rates"
