# #' Check data for functions
# #'
# #' @param X dfts or object that can be directly converted into dfts object
# #'
# #' @return dfts object or error if the data is incorrect
# #'
# #' @examples
# #' #.check_data(dfts(electricity))
# #' #.check_data(electricity)
# #'
# #' @keywords internal
# #' @noRd
# .check_data <- function(X, check.na=TRUE){
#   if(class(X)[1]=='dfts') {
#     if(check.na & sum(is.na(X$data))>0)
#       warning('NA values in data may affect some methods',call. = FALSE)
#     return(X)
#   }
#
#   if('data.frame' %in% class(X) || 'matrix' %in% class(X)){
#     if(check.na & sum(is.na(X))>0)
#       warning('NA values in data may affect some methods',call. = FALSE)
#     return(dfts(X))
#   }
#
#   stop('Check format of input data',call. = FALSE)
# }


#' Verify Inputs
#'
#' Standard function to do checks on input values
#'
#' @param selection User entered selection
#' @param poss_selections All possible sections
#'
#' @return The formatted input
#'
#' @keywords internal
#' @noRd
.verify_input <- function(selection, poss_selections){
  selection1 <- tolower(selection)
  poss_selections1 <- tolower(poss_selections)
  selection1 <- selection1[selection1 %in% poss_selections1]

  final_selection <- poss_selections[pmatch(selection1, poss_selections1)[1]]

  if(is.na(final_selection)){
    stop("Issue on the input '",deparse(substitute(selection)),
         "'. See documentation and verify your choice.",
         call. = FALSE)
  }

  final_selection
}
