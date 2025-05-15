#' plot_biplot
#'
#' Plot a PCA biplot of an assay in QFeatures object.
#'
#' @param qf
#' @param assay_name
#' @param color
#' @param shape
#' @param removeVar
#' @param legendPosition
#'
#' @returns
#' @export
#'
plot_biplot <- function(qf,
                        assay_name,
                        color = NULL,
                        shape = NULL,
                        removeVar = 0.1,
                        legendPosition = "right",
                        facet_formula = NULL) {

  # Extract assay and metadata
  se <- qf[[assay_name]]
  assay_data <- SummarizedExperiment::assay(se)
  meta <- as.data.frame(SummarizedExperiment::colData(se))

  # Dealing with weird shiny "" input problems
  if (is.null(color) || identical(color, "")) color <- NULL
  if (is.null(shape) || identical(shape, "")) shape <- NULL
  if (is.null(facet_formula) || identical(facet_formula, "")) facet_formula <- NULL

  # Run PCA
  p <- PCAtools::pca(assay_data,
                     metadata = meta,
                     removeVar = removeVar)

  # Combine PCA results and metadata
  df <- cbind(p$rotated, p$metadata)

  # Build dynamic aesthetics
  aes_args <- list(x = quote(PC1), y = quote(PC2))
  if (!is.null(color)) aes_args$color <- as.name(color)
  if (!is.null(shape)) aes_args$shape <- as.name(shape)

  # Build plot
  plt <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::labs(
      x = paste0("PC1 (", round(p$variance[1], 1), "%)"),
      y = paste0("PC2 (", round(p$variance[2], 1), "%)")
    ) +
    ggplot2::theme(legend.position = legendPosition)

  # Add facetting if requested
  if (!is.null(facet_formula)) {
    plt <- plt + ggplot2::facet_wrap(facet_formula)
  }

  return(plt)
}

