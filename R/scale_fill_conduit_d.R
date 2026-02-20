#' Add conduit colors to fill discrete values with ggplot2 scales
#'
#' @param ... Arguments passed on to \code{ggplot2::scale_fill_manual()}
#'
#' @returns ggproto object, scaleDiscrete
#' @export
#'
#' @examples
scale_fill_conduit_d <- function(...) {
  ggplot2::scale_fill_manual(
    values = rep(conduit_palette(), 100),
    ...
  )
}

