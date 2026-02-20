#' Plot Missing Value Pattern Heatmap
#'
#' Creates a heatmap visualization showing the pattern of missing values in a dataset,
#' with options for annotating both rows and columns using metadata variables.
#' This is useful for identifying systematic patterns in missing data and potential
#' batch effects or experimental issues.
#'
#' @param qf A QFeatures object containing the data to analyze
#' @param assay_name Character string specifying which assay to use for the analysis
#' @param col_color_variables Optional character vector specifying which columns from
#'   colData to use for column annotations (e.g., c("group", "batch"))
#' @param row_color_variables Optional character vector specifying which columns from
#'   rowData to use for row annotations (e.g., c("protein_class", "pathway"))
#' @param scale Logical indicating whether to scale the data before plotting
#'   (default: FALSE)
#'
#' @return A ComplexHeatmap object containing:
#'   \itemize{
#'     \item Main heatmap showing missing value pattern (black = present, white = missing)
#'     \item Column annotations using viridis color palette
#'     \item Row annotations using magma color palette
#'     \item Legend indicating missing and valid values
#'   }
#'
#' @export
#'
#' @examples
#' # Basic missing value heatmap:
#' # plot_missing_val_heatmap(qfeatures_obj, "protein")
#'
#' # With sample annotations:
#' # plot_missing_val_heatmap(qfeatures_obj, "protein",
#' #                         col_color_variables = c("group", "batch"))
#'
#' # With both sample and feature annotations:
#' # plot_missing_val_heatmap(qfeatures_obj, "protein",
#' #                         col_color_variables = c("treatment", "replicate"),
#' #                         row_color_variables = c("protein_class"))
plot_missing_val_heatmap <- function(qf,
                                     assay_name,
                                     col_color_variables = NULL,
                                     row_color_variables = NULL,
                                     max_rows = 5000) {

  # Extract SummarizedExperiment
  se <- qf[[assay_name]]
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # For taxonomic aggregations, treat 0 as NA
  assay_data <- assay(se)
  assay_data[assay_data == 0] <- NA

  ###### Missing value matrix ######
  missval <- +(!is.na(assay_data))  # 1 = present, 0 = missing

  # Safety check
  if (length(unique(missval)) < 2) {
    warning("All values are either missing or present. Heatmap cannot be plotted.")
    return(NULL)
  }

  ###### Downsample rows if needed ######
  n_rows <- nrow(missval)
  if (n_rows > max_rows) {
    message(glue::glue("Downsampling from {n_rows} to {max_rows} rows for visualization"))
    sampled_rows <- withr::with_seed(123, sample(n_rows, max_rows))

    # Subset the entire SE object, not just rowData
    se <- se[sampled_rows, , drop = FALSE]
    missval <- missval[sampled_rows, , drop = FALSE]
  }

  ###### Build annotation color maps dynamically ######
  metadata(se)$anno_colors <- list()

  if (!is.null(col_color_variables)) {
    for (var in col_color_variables) {
      values <- na.omit(unique(colData(se)[[var]]))
      metadata(se)$anno_colors[[var]] <- setNames(
        viridis::viridis(length(values)),
        values
      )
    }
  }

  if (!is.null(row_color_variables)) {
    for (var in row_color_variables) {
      values <- na.omit(unique(rowData(se)[[var]]))
      metadata(se)$anno_colors[[var]] <- setNames(
        viridis::magma(length(values)),
        values
      )
    }
  }

  ###### Create the heatmap ######
  ht <- ComplexHeatmap::Heatmap(
    missval,
    col = c("white", "black"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,  # improves performance on large matrices
    name = "Missing values",
    column_names_gp = grid::gpar(fontsize = 12),
    heatmap_legend_param = list(
      at = c(0, 1),
      labels = c("Missing", "Present")
    ),
    top_annotation = if (!is.null(col_color_variables))
      ComplexHeatmap::HeatmapAnnotation(
        df = colData(se)[, col_color_variables, drop = FALSE],
        col = metadata(se)$anno_colors
      ),
    right_annotation = if (!is.null(row_color_variables))
      ComplexHeatmap::rowAnnotation(
        df = rowData(se)[, row_color_variables, drop = FALSE],
        col = metadata(se)$anno_colors
      )
  )

  ###### Draw heatmap ######
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
}

