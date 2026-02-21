#' Conduit Color Palette
#'
#' Returns the official Conduit color palette (nine hex codes) used for logos and
#' consistent plotting (e.g. with `scale_color_conduit_d`, `scale_fill_conduit_d`).
#'
#' @return A character vector of 9 hex color codes.
#'
#' @export
#'
#' @examples
#' conduit_palette()
#' # Use in ggplot: scale_fill_manual(values = conduit_palette())
conduit_palette <- function(){
  cp <- c(
    "#F3B24A",  # orange
    "#31B2CC",  # blue
    "#4D4D4D",  # darkgray
    "#1F8A8C",  # teal
    "#FF6F61",  # coral
    "#A4C639",  # lime
    "#9B59B6",  # purple
    "#5DADE2",  # skyblue
    "#000000"   # black
  )
  return(cp)
}

