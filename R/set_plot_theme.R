#' Set Global Plot Theme
#'
#' Sets a global theme for all subsequent ggplot2 plots in the session.
#' This function provides a convenient way to maintain consistent plot styling
#' across multiple visualizations.
#'
#' @param theme_name Character string specifying the theme to use. Must be one of:
#'   \itemize{
#'     \item "ggprism": Clean, publication-ready theme from ggprism package
#'     \item "theme_classic": Classic theme with no background grid
#'     \item "theme_bw": Black and white theme
#'     \item "theme_dark": Dark theme
#'     \item "theme_void": Minimal theme with no axes or grid
#'     \item "theme_light": Light theme with gray background
#'     \item "theme_minimal": Minimal theme (default)
#'   }
#'
#' @return Invisibly returns NULL. The function sets the global theme for all
#'   subsequent ggplot2 plots in the session.
#'
#' @export
#'
#' @examples
#' # Default: minimal theme for all subsequent ggplots
#' set_plot_theme()
#'
#' \dontrun{
#' set_plot_theme("theme_classic")
#' set_plot_theme("ggprism")  # requires ggprism
#' }
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
