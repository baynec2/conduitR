#' Create an MA Plot from Limma Results
#'
#' Generates an MA plot to visualize differential expression results from limma
#' analysis. The plot shows average expression on the x-axis and log fold change
#' on the y-axis, making it easy to identify intensity-dependent fold change bias
#' that is not visible in a volcano plot.
#'
#' @param limma_stats A data frame containing limma analysis results, typically
#'   \code{perform_limma_analysis()$top_table}. Must contain columns
#'   \code{AveExpr}, \code{logFC}, and \code{adj.P.Val}.
#' @param color_by Optional character string specifying a column in
#'   \code{limma_stats} to use for point colors. When \code{NULL} (default),
#'   points are colored by significance (determined by \code{pval_threshold} and
#'   \code{logFC_threshold}).
#' @param pval_threshold Numeric significance threshold for adjusted p-values
#'   (default: 0.05). Used only when \code{color_by} is \code{NULL}.
#' @param logFC_threshold Numeric minimum absolute log fold change to consider
#'   significant (default: 1). Used only when \code{color_by} is \code{NULL}.
#' @param add_loess Logical; whether to overlay a loess smoothing curve to
#'   visualize intensity-dependent fold change bias (default: \code{TRUE}).
#' @param facet_formula Optional formula for faceting the plot (e.g., \code{~batch}).
#'
#' @return A ggplot object containing the MA plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' limma_res <- perform_limma_analysis(qf, "protein_groups", ~group, "groupB - groupA")
#' plot_ma(limma_res$top_table)
#' plot_ma(limma_res$top_table, pval_threshold = 0.01, logFC_threshold = 2)
#' plot_ma(limma_res$top_table, add_loess = FALSE)
#' plot_ma(limma_res$top_table, facet_formula = ~batch)
#' }
plot_ma <- function(limma_stats,
                    color_by = NULL,
                    pval_threshold = 0.05,
                    logFC_threshold = 1,
                    add_loess = TRUE,
                    facet_formula = NULL) {

  # NULL-guard for Shiny "" inputs
  if (is.null(color_by) || identical(color_by, "")) color_by <- NULL
  if (is.null(facet_formula) || identical(facet_formula, "")) facet_formula <- NULL

  if (!is.null(color_by)) {
    p1 <- limma_stats |>
      ggplot2::ggplot(ggplot2::aes(x = AveExpr, y = logFC,
                                   color = .data[[color_by]]))
  } else {
    limma_stats <- limma_stats |>
      dplyr::mutate(
        significant = adj.P.Val < pval_threshold & abs(logFC) > logFC_threshold
      )
    p1 <- limma_stats |>
      ggplot2::ggplot(ggplot2::aes(x = AveExpr, y = logFC, color = significant))
  }

  p1 <- p1 + ggplot2::geom_point(alpha = 0.5)

  # Reference line at logFC = 0
  p1 <- p1 + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                  color = "grey40")

  # Optional loess curve — inherit.aes = FALSE prevents per-group curves
  if (add_loess) {
    p1 <- p1 + ggplot2::geom_smooth(
      ggplot2::aes(x = AveExpr, y = logFC),
      method = "loess", se = FALSE,
      color = "steelblue", linewidth = 0.8,
      inherit.aes = FALSE
    )
  }

  if (!is.null(facet_formula)) {
    p1 <- p1 + ggplot2::facet_wrap(facet_formula)
  }

  p1 <- p1 + ggplot2::labs(x = "Average Expression", y = "Log Fold Change")

  return(p1)
}
