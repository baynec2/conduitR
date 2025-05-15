#' calc_relative_abundance
#'
#' calculate the relative abundance of an QFeatures assay.
#'
#' @param qf
#' @param assay_name
#'
#' @returns
#' @export
#'
#' @examples
calc_relative_abundance <- function(qf, assay_name) {

  # Check if the specified assay exists in the QFeatures object
  if (!(assay_name %in% names(qf))) {
    stop("The specified assay name does not exist in the QFeatures object.")
  }

  # Extract the assay data
  assay_data <- SummarizedExperiment::assay(qf, assay_name)

  # Calculate column sums (total abundance for each sample)
  col_sums <- colSums(assay_data, na.rm = TRUE)

  # Normalize each value by its corresponding column sum to get relative abundance
  rel_abundance <- sweep(assay_data, 2, col_sums, FUN = "/")

  # Ensure the relative abundance is wrapped in a SummarizedExperiment
  rel_abundance_se <- SummarizedExperiment::SummarizedExperiment(assays = list(rel_abundance =rel_abundance))

  # Add the relative abundance as a new assay in the QFeatures object
  qf <- QFeatures::addAssay(qf, rel_abundance_se, name = paste0(assay_name, "_rel_abundance"))

  # Return the modified QFeatures object
  return(qf)
}
