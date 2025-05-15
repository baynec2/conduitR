#' plot_sample_cor_heatmap
#'
#' plot the pearson correlation.
#'
#' @param qf
#' @param assay_name
#' @param sample_annotation_variables
#'
#' @returns
#' @export
#'
#' @examples
plot_sample_cor_heatmap = function(qf,
                                   assay_name,
                                   sample_annotation_variables){

  # Extracting relevant summarized experiment from QFeatues
  se = qf[[assay_name]]
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Generating the annotation colors
  metadata(se)$anno_colors <- list()
  shared_anno_colors <- list()

  # Generate consistent color palettes for each annotation variable
  for (var in sample_annotation_variables) {
    values <- unique(colData(se)[[var]])
    values <- values[!is.na(values)]
    palette <- viridis::viridis(length(values))
    shared_anno_colors[[var]] <- setNames(palette, values)
  }

  # Column annotations
  col_anno <- if (length(sample_annotation_variables) > 0) {
    ComplexHeatmap::HeatmapAnnotation(
      df = as.data.frame(colData(se)[, sample_annotation_variables, drop = FALSE]),
      col = shared_anno_colors
    )
  } else NULL

  # Row annotations
  row_anno <- if (length(sample_annotation_variables) > 0) {
    ComplexHeatmap::rowAnnotation(
      df = as.data.frame(colData(se)[, sample_annotation_variables, drop = FALSE]),
      col = shared_anno_colors
    )
  } else NULL

  # Generating a correlation matrix
  cor_mat <- cor(assay(se), method = "pearson")

  # Setting the color palette
  col = viridis::viridis(1000)

  # Plotting
  p1 = ComplexHeatmap::Heatmap(
    cor_mat,
    col = col,
    name = "Correlation",
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = col_anno,
    right_annotation = row_anno
  )

  return(p1)
}
