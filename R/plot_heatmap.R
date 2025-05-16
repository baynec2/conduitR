#' Create Static Heatmap with Annotations
#'
#' Generates a static heatmap visualization using the sechm package,
#' with options for row and column annotations and hierarchical clustering.
#' This function provides a simpler alternative to plot_heatmaply for
#' non-interactive visualizations.
#'
#' @param qf A QFeatures object containing the data to visualize
#' @param assay_name Character string specifying which assay to use for the heatmap
#' @param col_color_variables Optional character vector specifying which columns from
#'   colData to use for column annotations (e.g., c("group", "batch"))
#' @param row_color_variables Optional character vector specifying which columns from
#'   rowData to use for row annotations (e.g., c("protein_class", "pathway"))
#' @param scale Logical indicating whether to scale the data (default: TRUE)
#' @param ... Additional arguments passed to sechm::sechm()
#'
#' @return A ComplexHeatmap object containing:
#'   \itemize{
#'     \item Main heatmap showing the data matrix
#'     \item Column annotations using viridis color palette
#'     \item Row annotations using magma color palette
#'     \item Hierarchical clustering of both rows and columns
#'     \item Color scale from viridis palette
#'   }
#'
#' @export
#'
#' @examples
#' # Basic heatmap:
#' # plot_heatmap(qfeatures_obj, "protein")
#' 
#' # With sample annotations:
#' # plot_heatmap(qfeatures_obj, "protein",
#' #             col_color_variables = c("group", "batch"))
#' 
#' # With both sample and feature annotations:
#' # plot_heatmap(qfeatures_obj, "protein",
#' #             col_color_variables = c("treatment", "replicate"),
#' #             row_color_variables = c("protein_class"))
#' 
#' # Without scaling:
#' # plot_heatmap(qfeatures_obj, "protein", scale = FALSE)
plot_heatmap = function(qf,
                        assay_name,
                        col_color_variables = NULL,
                        row_color_variables = NULL,
                        scale = TRUE,
                        ...) {

  # Load required package
  requireNamespace("viridis", quietly = TRUE)

  # Extracting assay from QF object
  se <- qf[[assay_name]]

  # Set heatmap colors
  metadata(se)$hmcols <- viridis::viridis(1000)

  # Initialize color list
  metadata(se)$anno_colors <- list()

  # Assign column annotation colors dynamically
  if (!is.null(col_color_variables)) {
    for (var in col_color_variables) {
      values <- unique(colData(se)[[var]])
      values <- values[!is.na(values)]
      n_values <- length(values)
      palette <- viridis::viridis(n_values)
      metadata(se)$anno_colors[[var]] <- setNames(palette, values)
    }
  }

  # Assign row annotation colors dynamically
  if (!is.null(row_color_variables)) {
    for (var in row_color_variables) {
      values <- unique(rowData(se)[[var]])
      values <- values[!is.na(values)]
      n_values <- length(values)
      palette <- viridis::magma(n_values)
      metadata(se)$anno_colors[[var]] <- setNames(palette, values)
    }
  }

  # Generate heatmap
  p1 <- sechm::sechm(
    se,
    features = rownames(se),
    top_annotation = col_color_variables,
    left_annotation = row_color_variables,
    do.scale = scale,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    ...
  )

  return(p1)
}

