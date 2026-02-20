#' Conduit color palette
#'
#' @returns A character vector of hex color codes representing the Conduit color
#' palette.
#' @export
#'
#' @examples
#' # Returns 9 colors that go along with the conduit logos.
#' conduit_palette()
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

