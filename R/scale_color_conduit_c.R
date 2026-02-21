#' Add Conduit Colors for Continuous Color in ggplot2
#'
#' Wraps `ggplot2::scale_color_gradient()` with conduit endpoints (blue low,
#' orange high). Use for continuous color aesthetics (e.g. point size/color).
#'
#' @param ... Arguments passed on to `ggplot2::scale_color_gradient()`.
#'
#' @return A `ScaleContinuous` (ggproto) object for color.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = mpg, y = wt, color = disp)) +
#'   geom_point() +
#'   scale_color_conduit_c()
#' }
scale_color_conduit_c <- function(...) {
  ggplot2::scale_color_gradient(
    low  = "#31B2CC",  # blue = low values
    high = "#F3B24A",  # orange = high values
    ...
  )
}
