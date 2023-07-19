#' Spanish Spot Electricity Data
#'
#' A fChange ready data set containing the hourly electricity spot prices from
#'  Spain in 2014.
#'
#' @format ## `electricity`
#' A data frame with 24 rows and 365 columns, representing the hourly (row)
#'  electrity spot values for each day (column) of the year in 2014.
#'
#'@source <www.omie.es>
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
#'@source <www.bom.gov.au>
"temperature"
