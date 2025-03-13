#' prepare_qfeature
#'
#' assemble qfeature object from a list of any number of matrices corresponding
#' to the outputs of the process.
#'
#' @param sample_annotation_fp filepath of sample annotation file path
#' @param vector_of_matrix_fps vector of filepath
#' @param vector_of_matrices_names vector of names to call each matrices.
#' Optional, by default the file name - the extension will be used.
#' @returns
#' @export
#'
#' @examples
#'
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
