#' Prepare QFeatures Object from Multiple Data Matrices
#'
#' Assembles a QFeatures object from a list of data matrices and sample annotations.
#' This function is particularly useful for integrating multiple types of quantitative
#' data (e.g., protein abundances, peptide intensities) into a single, organized
#' data structure that maintains relationships between features and samples.
#'
#' @param sample_annotation_fp Character string specifying the file path to the sample
#'   annotation file. This file should be a tab-delimited text file containing sample
#'   metadata, with a "file" column that matches the column names in the data matrices.
#' @param vector_of_matrix_fps Character vector containing file paths to the data matrices.
#'   Each file should be a tab-delimited text file with features as rows and samples
#'   as columns. The first column should contain feature identifiers.
#' @param vector_of_matrices_names Optional character vector specifying names for each
#'   matrix in the QFeatures object. If NULL (default), names are derived from the
#'   matrix filenames by removing "_matrix.tsv" extension.
#'
#' @return A QFeatures object containing:
#'   \itemize{
#'     \item Multiple assays, one for each input matrix
#'     \item Shared sample annotations (colData) across all assays
#'     \item Feature annotations (rowData) specific to each assay
#'     \item Consistent sample identifiers across all components
#'   }
#'   Each assay is stored as a SummarizedExperiment object with:
#'   \itemize{
#'     \item assay: The quantitative data matrix
#'     \item rowData: Feature metadata from the first column of each matrix
#'     \item colData: Sample metadata from the annotation file
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage with default matrix names:
#' # qf <- prepare_qfeature(
#' #   sample_annotation_fp = "samples.tsv",
#' #   vector_of_matrix_fps = c("protein_matrix.tsv", "peptide_matrix.tsv")
#' # )
#' 
#' # Specify custom names for the matrices:
#' # qf <- prepare_qfeature(
#' #   sample_annotation_fp = "samples.tsv",
#' #   vector_of_matrix_fps = c("protein_matrix.tsv", "peptide_matrix.tsv"),
#' #   vector_of_matrices_names = c("protein", "peptide")
#' # )
#' 
#' # The resulting QFeatures object can be used with other functions:
#' # plot_heatmap(qf, "protein")
#' # calc_relative_abundance(qf, "protein")
#' # perform_limma_analysis(qf, "protein", ~group, "groupB - groupA")
#'
#' @note
#' This function:
#' \itemize{
#'   \item Requires consistent sample identifiers across all inputs
#'   \item Automatically handles feature identifiers as row names
#'   \item Preserves all sample metadata from the annotation file
#'   \item Creates a hierarchical data structure for multi-level analysis
#' }
#' 
#' Important considerations:
#' \itemize{
#'   \item Sample annotation file must have a "file" column
#'   \item Matrix files must be tab-delimited
#'   \item Column names in matrices must match sample identifiers
#'   \item Feature IDs should be unique within each matrix
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[QFeatures]{QFeatures}} for the underlying data structure
#'   \item \code{\link[SummarizedExperiment]{SummarizedExperiment}} for assay structure
#'   \item \code{\link[plot_heatmap]{plot_heatmap}} for visualizing the data
#'   \item \code{\link[calc_relative_abundance]{calc_relative_abundance}} for
#'     calculating relative abundances
#'   \item \code{\link[perform_limma_analysis]{perform_limma_analysis}} for
#'     differential analysis
#' }
prepare_qfeature <- function(sample_annotation_fp,
                             vector_of_matrix_fps,
                             vector_of_matrices_names = NULL) {

  # Read the annotation file
  colData <- readr::read_delim(sample_annotation_fp) |>
    tibble::column_to_rownames("file")  # Use file as sample identifier

  quantCols <- rownames(colData)  # Sample identifiers

  # Load matrices into list
  matrix_list <- purrr::map(vector_of_matrix_fps, readr::read_tsv)

  # Assign names to matrices (if not provided, derive from filenames)
  names(matrix_list) <- if (is.null(vector_of_matrices_names)) {
    gsub("_matrix.tsv", "", basename(vector_of_matrix_fps))
  } else {
    vector_of_matrices_names
  }

  # Convert each matrix into a SummarizedExperiment object
  se_list <- purrr::map(matrix_list, function(mat) {

    # Convert first column to row names (features), no need to remove it
    rowData <- mat |>
      tibble::column_to_rownames(colnames(mat)[1]) |>
      dplyr::select(-quantCols)# Convert first column (feature IDs) to row names

    # Extract the numeric data (assay matrix), excluding the first column (already handled)
    assay_matrix <- mat |>
      tibble::column_to_rownames(colnames(mat)[1]) |>
      dplyr::select(quantCols)# Numeric data (no feature IDs)

    # Ensure colnames match sample_annotation rownames
    if (!all(rownames(colData) %in% colnames(assay_matrix))) {
      stop("Mismatch: Column names in matrix do not match sample annotation.")
    }

    # Create SummarizedExperiment object
    SummarizedExperiment::SummarizedExperiment(
      assays = list(intensity = assay_matrix),  # Assay data (counts)
      colData = colData,  # Match colData
      rowData = rowData  # Match rowData
    )
  })

  # Generate QFeatures object
  QF <- QFeatures::QFeatures(se_list,colData = colData)

  return(QF)
}
