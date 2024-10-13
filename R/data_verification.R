#' Check data for functions
#'
#' @param X funts or object that can be directly converted into funts object
#'
#' @return funts object or error if the data is incorrect
#'
#' @examples
#' .check_data(funts(electricity))
#' .check_data(electricity)
#'
#' @keywords internal
.check_data <- function(X, check.na=TRUE){
  if(class(X)[1]=='funts') {
    if(check.na & sum(is.na(X$data))>0)
      warning('NA values in data may affect some methods',call. = FALSE)
    return(X)
  }

  if('data.frame' %in% class(X) || 'matrix' %in% class(X)){
    if(check.na & sum(is.na(X))>0)
      warning('NA values in data may affect some methods',call. = FALSE)
    return(funts(X))
  }

  stop('Check format of input data',call. = FALSE)
}


#' Title
#'
#' @param selection
#' @param poss_selections
#'
#' @return
#' @export
#'
#' @examples
#'
#' @keywords internal
.verify_input <- function(selection, poss_selections){
  # TODO:: Use this in functions instead of repeating code
  final_selection <- poss_selections[min(pmatch(selection, poss_selections))]

  if(is.na(final_selection)){
    stop("Issue on the input '",deparse(substitute(selection)),
         "'. See documentation and verify your choice.",
         call. = FALSE)
  }
  final_selection
}
