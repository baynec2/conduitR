#' Negation of `%in%`
#'
#' This operator returns `TRUE` for elements not in the second argument.
#'
#' @param x Vector or NULL: the values to be matched.
#' @param table Vector or NULL: the values to be matched against.
#' @return A logical vector.
#' @export
`%!in%` <- Negate(`%in%`)
