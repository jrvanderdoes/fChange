#' Forecasting competition Evaluation
#'
#' @param estimate Matrix of forecast
#'
#' @returns MSE for competition
#'
#' @noRd
#' @keywords internal
competition_error <- function(estimate){

  urlfile <-  paste0('https://raw.githubusercontent.com/jrvanderdoes/fChange',
                     '/main/vignettes/articles/functionaldata/practice.rds')
  data_true <- readRDS(urlfile)

  if(is.dfts(estimate)) estimate <- estimate$data

  error <- rep(NA, 14)
  for(i in 1:14){
    error[i] <- sum((estimate[,i] - data_true$data[,i])^2)
  }

  mean(error)
}
