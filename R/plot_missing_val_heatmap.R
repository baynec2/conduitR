#' plot_missing_val_heatmap
#'
#' @param qf
#' @param assay_name
#' @param col_color_variables
#' @param row_color_variables
#' @param scale
#'
#' @returns
#' @export
#'
#' @examples
plot_missing_val_heatmap = function(qf,
                                    assay_name,
                                    col_color_variables = NULL,
                                    row_color_variables = NULL,
                                    scale = FALSE) {

  # Extracting the SummarizedExperiment
  se = qf[[assay_name]]

  # Ensure the object is a SummarizedExperiment
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Allow the user the ability to annotate the rows and columns by colors
  metadata(se)$anno_colors <- list()

  # Assign column annotation colors dynamically if specified
  col_anno <- NULL
  if (!is.null(col_color_variables)) {
    for (var in col_color_variables) {
      values <- unique(colData(se)[[var]])
      values <- values[!is.na(values)]  # Removing NA values
      n_values <- length(values)
      palette <- viridis::viridis(n_values)
      metadata(se)$anno_colors[[var]] <- setNames(palette, values)

      # Prepare for ComplexHeatmap column annotation
      col_anno <- c(col_anno, var)
    }
  }

  # Assign row annotation colors dynamically if specified
  row_anno <- NULL
  if (!is.null(row_color_variables)) {
    for (var in row_color_variables) {
      values <- unique(rowData(se)[[var]])
      values <- values[!is.na(values)]  # Removing NA values
      n_values <- length(values)
      palette <- viridis::magma(n_values)
      metadata(se)$anno_colors[[var]] <- setNames(palette, values)

      # Prepare for ComplexHeatmap row annotation
      row_anno <- c(row_anno, var)
    }
  }

  # Extract the assay data
  se_assay <- assay(se)

  # Find rows with missing values
  missval <- ifelse(is.na(se_assay), 0, 1)

  # Create the heatmap with missing value pattern
  ht2 = ComplexHeatmap::Heatmap(
    missval,
    col = c("white", "black"),
    column_names_side = "top",
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "Missing values pattern",
    column_names_gp = grid::gpar(fontsize = 16),
    heatmap_legend_param = list(
      at = c(0, 1),
      labels = c("Missing value", "Valid value")
    ),
    top_annotation = if (!is.null(col_anno)) ComplexHeatmap::HeatmapAnnotation(col_annot = colData(se)[, col_anno]),
    left_annotation = if (!is.null(row_anno)) ComplexHeatmap::HeatmapAnnotation(row_annot = rowData(se)[, row_anno])
  )

  # Draw the heatmap with the legend on top
  ComplexHeatmap::draw(ht2, heatmap_legend_side = "top")
}
