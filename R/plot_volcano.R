#' plot_volcano()
#'
#'
#'
#' @param limma_stats
#' @param facet_formula
#' @param color_by
#' @param pval_threshold
#'
#' @returns
#' @export
#'
#' @examples
plot_volcano <- function(limma_stats,
                         facet_formula = NULL,
                         color_by = NULL,
                         pval_threshold = 0.05) {

  p1 <- limma_stats |>
    ggplot2::ggplot(ggplot2::aes(x = logFC, y = neg_log10.adj.P.Value))

  # Dealing with weird shiny "" input problems
  if (is.null(color_by) || identical(color_by, "")) color_by <- NULL
  if (is.null(facet_formula) || identical(facet_formula, "")) facet_formula <- NULL

  # Add colored points if specified
  if (!is.null(color_by)) {
    p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]))
  } else {
    p1 <- p1 + ggplot2::geom_point()
  }

  # Facet if formula is given
  if (!is.null(facet_formula) && facet_formula != "") {
    p1 <- p1 + ggplot2::facet_wrap(facet_formula)
  }

  # Add horizontal line for significance threshold
  p1 <- p1 + ggplot2::geom_hline(yintercept = -log10(pval_threshold),
                                 linetype = "dashed",
                                 color = "red")

  return(p1)
}


