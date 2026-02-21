#' Add Conduit Colors for Discrete Color in ggplot2
#'
#' Wraps `ggplot2::scale_color_manual()` with the conduit color palette (recycled
#' for many levels). Use for coloring points/lines by a discrete variable.
#'
#' @param ... Arguments passed on to `ggplot2::scale_color_manual()` (e.g.
#'   `name`, `na.value`).
#'
#' @return A `ScaleDiscrete` (ggproto) object for color.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = mpg, y = wt, color = factor(cyl))) +
#'   geom_point() +
#'   scale_color_conduit_d()
#' }
scale_color_conduit_d <- function(...) {
  ggplot2::scale_color_manual(
    values = rep(conduit_palette(), 100),  # recycle colors for any number of levels
    ...
  )
}
