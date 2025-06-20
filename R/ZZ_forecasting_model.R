#' Forecasting competition Evaluation
#'
#' @param estimate Matrix of forecast (matrix of 24x30)
#'
#' @returns MSE for competition
#'
#' @noRd
#' @keywords internal
competition_error <- function(estimate){

  urlfile <-  paste0('https://raw.githubusercontent.com/jrvanderdoes/fChange',
                     '/main/vignettes/articles/functionaldata/practice.csv')
  data_true <- utils::read.csv(urlfile)[,-1]

  if(is.dfts(estimate)) estimate <- estimate$data

  error <- rep(NA, 14)
  for(i in 1:14){
    error[i] <- sum((estimate[,i] - data_true[,i])^2)
  }

  mean(error)
}
