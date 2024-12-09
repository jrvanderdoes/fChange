#' Pipe operator
#'
#' See \code{tidyr::\link[tidyr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#'
#' @rdname pipe
#'
#' @keywords internal
#' @export
#'
#' @importFrom tidyr %>%
#'
#' @usage lhs \%>\% rhs
#'
#' @param lhs A value or the tidyr placeholder.
#' @param rhs A function call using the tidyr semantics.
#'
#' @return The result of calling \code{rhs(lhs)}.
NULL
