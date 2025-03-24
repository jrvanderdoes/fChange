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
  # method <- match.arg(method, c('welch','mc','imhof'))

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
