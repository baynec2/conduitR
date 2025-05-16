#' Plot Sample Correlation Heatmap
#'
#' Creates a heatmap visualization of Pearson correlations between samples,
#' with optional sample annotations displayed as colored bars. This is useful
#' for assessing sample relationships and identifying potential batch effects
#' or experimental groups.
#'
#' @param qf A QFeatures object containing the data to analyze
#' @param assay_name Character string specifying which assay to use for correlation analysis
#' @param sample_annotation_variables Character vector specifying which columns from
#'   colData to use for sample annotations (e.g., c("group", "batch", "treatment"))
#'
#' @return A ComplexHeatmap object containing:
#'   \itemize{
#'     \item Main heatmap showing Pearson correlations between samples
#'     \item Color scale from viridis palette
#'     \item Top and right annotations showing sample metadata
#'     \item Consistent color schemes for annotation variables
#'   }
#'
#' @export
#'
#' @examples
#' # Basic correlation heatmap:
#' # plot_sample_cor_heatmap(qfeatures_obj, "protein", NULL)
#' 
#' # With sample annotations:
#' # plot_sample_cor_heatmap(qfeatures_obj, "protein", 
#' #                        c("group", "batch"))
#' 
#' # With multiple annotation variables:
#' # plot_sample_cor_heatmap(qfeatures_obj, "protein",
#' #                        c("treatment", "timepoint", "replicate"))
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
