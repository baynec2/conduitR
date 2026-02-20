#'  Add conduit colors to color discrete values with ggplot2 scales
#'
#' @param ... Arguments passed on to \code{ggplot2::scale_color_manual()}
#'
#' @returns ggproto object, scaleDiscrete
#' @export
#'
#' @examples
scale_color_conduit_d <- function(...) {
  ggplot2::scale_color_manual(
    values = rep(conduit_palette(), 100),  # recycle colors for any number of levels
    ...
  )
}
