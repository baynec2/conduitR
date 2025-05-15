#' plot_heatmap
#'
#' Easy heatmap plots.
#'
#' @param qf
#' @param assay_name
#' @param col_color_variables
#' @param row_color_variables
#' @param scale
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
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

