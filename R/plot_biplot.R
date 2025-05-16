#' Create PCA Biplot from QFeatures Data
#'
#' Generates a PCA biplot visualization for a specified assay in a QFeatures object,
#' allowing for customization of point aesthetics and faceting. This function uses
#' the PCAtools package for PCA computation and visualization.
#'
#' @param qf A QFeatures object containing the assay data to plot
#' @param assay_name Character string specifying which assay to use for the PCA
#' @param color Optional character string specifying a column in the metadata to use for point colors
#' @param shape Optional character string specifying a column in the metadata to use for point shapes
#' @param removeVar Numeric value between 0 and 1 specifying the proportion of variance to remove
#'   before PCA (default: 0.1)
#' @param legendPosition Character string specifying legend position (default: "right")
#' @param facet_formula Optional formula for faceting the plot (e.g., ~group)
#'
#' @return A ggplot object containing the PCA biplot with:
#'   \itemize{
#'     \item Points colored and/or shaped by specified metadata variables
#'     \item Axis labels showing percentage of variance explained
#'     \item Optional faceting based on the provided formula
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage:
#' # plot_biplot(qfeatures_obj, "protein")
#' 
#' # With aesthetics:
#' # plot_biplot(qfeatures_obj, "protein", color = "group", shape = "treatment")
#' 
#' # With faceting:
#' # plot_biplot(qfeatures_obj, "protein", facet_formula = ~group)
#'
#' @note
#' This function requires the PCAtools package. For large datasets (e.g., >10,000 features),
#' consider using removeVar to reduce dimensionality before PCA computation to improve
#' performance and memory usage.
#'
#' @seealso \code{\link[PCAtools]{pca}} for the underlying PCA computation
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

