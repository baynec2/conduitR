#' Create a Volcano Plot from Limma Results
#'
#' Generates a volcano plot to visualize differential expression results from limma analysis.
#' The plot shows log fold changes on the x-axis and -log10 adjusted p-values on the y-axis,
#' with options for coloring points and faceting.
#'
#' @param limma_stats A data frame containing limma analysis results, typically the output
#'   from perform_limma_analysis(). Must contain columns 'logFC' and 'neg_log10.adj.P.Value'.
#' @param facet_formula Optional formula for faceting the plot (e.g., ~group)
#' @param color_by Optional character string specifying a column to use for point colors
#' @param pval_threshold Numeric value specifying the significance threshold for p-values
#'   (default: 0.05). A horizontal line will be drawn at -log10(pval_threshold).
#'
#' @return A ggplot object containing the volcano plot with:
#'   \itemize{
#'     \item Points representing proteins/features
#'     \item X-axis showing log fold change
#'     \item Y-axis showing -log10 adjusted p-value
#'     \item Optional coloring by specified variable
#'     \item Optional faceting
#'     \item Horizontal line indicating significance threshold
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' limma_res <- perform_limma_analysis(qf, "protein_groups", ~ group, "B - A")
#' plot_volcano(limma_res$top_table)
#' plot_volcano(limma_res$top_table, color_by = "group", pval_threshold = 0.01)
#' plot_volcano(limma_res$top_table, facet_formula = ~ batch)
#' }
plot_volcano <- function(limma_stats,
                         facet_formula = NULL,
                         color_by = NULL,
                         pval_threshold = 0.05) {

  p1 <- limma_stats |>
    ggplot2::ggplot(ggplot2::aes(x = logFC, y = neg_log10.adj.P.Val))

  # Dealing with weird shiny "" input problems
  if (is.null(color_by) || identical(color_by, "")) color_by <- NULL
  if (is.null(facet_formula) || identical(facet_formula, "")) facet_formula <- NULL

  # Add colored points if specified
  if (!is.null(color_by)) {
    p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]),alpha = 0.5)
  } else {
    p1 <- p1 + ggplot2::geom_point(alpha = 0.5)
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


