#' Create Interactive Heatmap with Annotations
#'
#' Generates an interactive heatmap visualization using the heatmaply package,
#' allowing for interactive exploration of data with hover information and
#' optional row and column annotations.
#'
#' @param qf A QFeatures object containing the data to visualize
#' @param assay_name Character string specifying which assay to use for the heatmap
#' @param col_color_variables Optional character vector specifying which columns from
#'   colData to use for column annotations (e.g., c("group", "batch"))
#' @param row_color_variables Optional character vector specifying which columns from
#'   rowData to use for row annotations (e.g., c("protein_class", "pathway"))
#' @param scale Character string specifying how to scale the data. Must be one of:
#'   \itemize{
#'     \item "column": Scale each column (default)
#'     \item "row": Scale each row
#'     \item "none": No scaling
#'   }
#' @param hover_col_variables Optional character vector specifying which columns from
#'   colData to show in hover text. Use "all" to show all columns.
#' @param hover_row_variables Optional character vector specifying which columns from
#'   rowData to show in hover text. Use "all" to show all columns.
#' @param ... Additional arguments passed to heatmaply::heatmaply()
#'
#' @return An interactive heatmaply object containing:
#'   \itemize{
#'     \item Main heatmap showing the data matrix
#'     \item Column annotations using specified variables
#'     \item Row annotations using specified variables
#'     \item Interactive hover information showing metadata
#'     \item Zoom and pan capabilities
#'   }
#'
#' @export
#'
#' @examples
#' # Basic interactive heatmap:
#' # plot_heatmaply(qfeatures_obj, "protein")
#' 
#' # With sample annotations and hover info:
#' # plot_heatmaply(qfeatures_obj, "protein",
#' #               col_color_variables = c("group", "batch"),
#' #               hover_col_variables = "all")
#' 
#' # With both sample and feature annotations:
#' # plot_heatmaply(qfeatures_obj, "protein",
#' #               col_color_variables = c("treatment", "replicate"),
#' #               row_color_variables = c("protein_class"),
#' #               hover_row_variables = c("description", "pathway"))
plot_heatmaply <- function(qf,
                           assay_name,
                           col_color_variables = NULL,
                           row_color_variables = NULL,
                           scale = "column",
                           hover_col_variables = NULL,
                           hover_row_variables = NULL,
                           ...) {

  # Extract SummarizedExperiment
  se <- qf[[assay_name]]

  # Extract expression matrix
  mat <- as.matrix(SummarizedExperiment::assay(se))

  # Column annotation (as factors)
  column_side_colors <- NULL
  if (!is.null(col_color_variables)) {
    column_side_colors <- as.data.frame(SummarizedExperiment::colData(se)[, col_color_variables, drop = FALSE])
    column_side_colors[] <- lapply(column_side_colors, as.factor)
  }

  # Row annotation (as factors)
  row_side_colors <- NULL
  if (!is.null(row_color_variables)) {
    row_side_colors <- as.data.frame(SummarizedExperiment::rowData(se)[, row_color_variables, drop = FALSE])
    row_side_colors[] <- lapply(row_side_colors, as.factor)
  }

  # 🔍 Custom hover text using colData and rowData
  coldata_full <- as.data.frame(SummarizedExperiment::colData(se))
  rowdata_full <- as.data.frame(SummarizedExperiment::rowData(se))

  if (is.null(hover_col_variables)) {
    coldata_hover <- NULL
  } else if (identical(hover_col_variables, "all")) {
    coldata_hover <- coldata_full
  } else if (all(hover_col_variables %in% colnames(coldata_full))) {
    coldata_hover <- coldata_full[, hover_col_variables, drop = FALSE]
  } else {
    stop("hover_col_variables must be NULL, 'all', or a vector of valid column names from colData.")
  }

  if (is.null(hover_row_variables)) {
    rowdata_hover <- NULL
  } else if (identical(hover_row_variables, "all")) {
    rowdata_hover <- rowdata_full
  } else if (all(hover_row_variables %in% colnames(rowdata_full))) {
    rowdata_hover <- rowdata_full[, hover_row_variables, drop = FALSE]
  } else {
    stop("hover_row_variables must be NULL, 'all', or a vector of valid column names from rowData.")
  }

  hover_text <- outer(
    rownames(mat),
    colnames(mat),
    Vectorize(function(row, col) {
      row_info <- paste(paste(names(rowdata_hover), rowdata_hover[row, ], sep = ": "), collapse = "<br>")
      col_info <- paste(paste(names(coldata_hover), coldata_hover[col, ], sep = ": "), collapse = "<br>")
      paste0("<br><br><b>Row Metadata</b><br>", row_info, "<br><br><b>Column Metadata</b><br>", col_info)
    })
  )

  # Plot
  p <- heatmaply::heatmaply(
    mat,
    scale = scale,
    showticklabels = c(FALSE, FALSE),
    col_side_colors = column_side_colors,
    row_side_colors = row_side_colors,
    custom_hovertext = hover_text,
    ...
  )

  return(p)
}
