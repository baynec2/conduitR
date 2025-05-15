#' set_plot_theme
#'
#' globally set the plot theme
#'
#' @param theme_name
#'
#' @returns
#' @export
#'
#' @examples
set_plot_theme <- function(theme_name = "theme_minimal") {
  `%!in%` <- Negate(`%in%`)

  supported_themes <- c(
    "ggprism",
    "theme_classic",
    "theme_bw",
    "theme_dark",
    "theme_void",
    "theme_light",
    "theme_minimal"
  )

  if (theme_name %!in% supported_themes) {
    stop("Error: Select a supported plot theme")
  }

  if (theme_name == "ggprism") {
    ggplot2::theme_set(ggprism::theme_prism())
  } else {
    ggplot2::theme_set(get(theme_name, envir = asNamespace("ggplot2"))())
  }
}
