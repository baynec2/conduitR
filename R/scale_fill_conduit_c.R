#' Add Conduit Colors for Continuous Fill in ggplot2
#'
#' Wraps `ggplot2::scale_fill_gradient()` with conduit palette endpoints (blue
#' low, orange high). Use with continuous fill aesthetics (e.g. heatmaps).
#'
#' @param ... Arguments passed on to `ggplot2::scale_fill_gradient()` (e.g.
#'   `name`, `na.value`).
#'
#' @return A `ScaleContinuous` (ggproto) object for fill.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = disp)) +
#'   geom_tile() +
#'   scale_fill_conduit_c()
#' }
scale_fill_conduit_c <- function(...) {
  ggplot2::scale_fill_gradient(
    low  = "#31B2CC",  # blue = low values
    high = "#F3B24A",  # orange = high values
    ...
  )
}
