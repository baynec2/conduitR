#' Add Conduit Colors for Discrete Fill in ggplot2
#'
#' Wraps `ggplot2::scale_fill_manual()` with the conduit palette (recycled for
#' many levels). Use for filling bars/areas by a discrete variable.
#'
#' @param ... Arguments passed on to `ggplot2::scale_fill_manual()`.
#'
#' @return A `ScaleDiscrete` (ggproto) object for fill.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = factor(cyl), fill = factor(gear))) +
#'   geom_bar() +
#'   scale_fill_conduit_d()
#' }
scale_fill_conduit_d <- function(...) {
  ggplot2::scale_fill_manual(
    values = rep(conduit_palette(), 100),
    ...
  )
}

