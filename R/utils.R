#' Negation of the `%in%` Operator
#'
#' A utility operator that returns the logical negation of the `%in%` operator.
#' This is useful for checking if elements are NOT present in a vector or set.
#'
#' @param x Vector or NULL: the values to be matched
#' @param table Vector or NULL: the values to be matched against
#'
#' @return A logical vector of the same length as x, indicating for each element
#'   whether it is NOT present in table. TRUE indicates the element is not in
#'   table, FALSE indicates it is present.
#'
#' @export
#'
#' @examples
#' # Check if elements are not in a vector
#' c("a", "b", "c") %!in% c("a", "d")
#' # Returns: FALSE TRUE TRUE
#'
#' # Filter data frame rows
#' df <- data.frame(id = c("a", "b", "c"), val = 1:3)
#' df[df$id %!in% c("a"), ]
#'
#' # Remove specific elements from a vector
#' x <- c("keep", "remove1", "keep2", "remove2")
#' x[x %!in% c("remove1", "remove2")]
#'
#' @seealso \code{\link[base]{\%in\%}} for the original operator
`%!in%` <- Negate(`%in%`)
