#'  Add conduit colors to color continuous values with ggplot2 scales
#'
#' @param ... Arguments passed on to \code{ggplot2::scale_color_gradient()}
#'
#' @returns ggproto object, scaleContinuous
#' @export
#'
#' @examples
scale_color_conduit_c <- function(...) {
  ggplot2::scale_color_gradient(
    low  = "#31B2CC",  # blue = low values
    high = "#F3B24A",  # orange = high values
    ...
  )
}
