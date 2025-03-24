#' Spanish Spot Electricity Data
#'
#' The hourly electricity spot prices from Spain in 2014.
#'
#' @format ## `electricity`
#' A dfts object.
#'
#' @source <www.omie.es>
"electricity"

#' Australian Temperature Data
#'
#' A list of daily minimum temperature data (1859 - 2012) from reporting station
#'  in Sydney, Australia. Any missing data was previously, linearly,
#'  interpolated and leap days were removed for consistency.
#'
#' @format ## `temperature`
#' A dfts object.
#'
#' @source <www.bom.gov.au>
"temperature"

#' S&P 500 Index Data
#'
#' Intraday prices for the S&P500 index (SPY) for 2019 to 2023 with holidays and
#'  weekends removed. Minutely resolution and daily observations.
#'
#' @format ## `SPY500`
#' A dfts object.
"SPYUS500"


#' Breast Cancer
#'
#' Percentage of cause-specific deaths out of total deaths for female breast
#'  cancer in the United States from 1950 to 2021.
#'
#' @format ## `cancer`
#' A data.frame with columns being the year and rows the age groups (5 years).
"cancer"


#' US Yield Curves
#'
#' Yield curves in the US for 1 - 360 month maturity from 1990 to 2023
#'  (removing days without information, i.e. weekends and holidays).
#'
#' @format ## `rates`
#' A dfts object.
"rates"
