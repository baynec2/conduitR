#' Calculate Relative Abundance for a QFeatures Assay
#'
#' Calculates the relative abundance of features (e.g., proteins or species) within each
#' sample of a QFeatures assay by normalizing the raw abundance values to proportions.
#' This function adds a new assay containing the relative abundance values to the
#' QFeatures object, making it easy to compare the composition of features across samples.
#'
#' @param qf A QFeatures object containing the abundance data. The assay should contain
#'   non-negative abundance values (e.g., protein intensities or species counts).
#' @param assay_name Character string specifying which assay to use for calculating
#'   relative abundances. The assay must exist in the QFeatures object.
#'
#' @return A QFeatures object with an additional assay named "{assay_name}_rel_abundance"
#'   containing:
#'   \itemize{
#'     \item Values normalized to proportions (0-1)
#'     \item Same dimensions as the input assay
#'     \item NA values preserved from the input
#'     \item Column sums equal to 1 for each sample (excluding NA values)
#'   }
#'   The original assays remain unchanged.
#'
#' @export
#'
#' @examples
#' # Calculate relative abundance for a protein assay:
#' # qf_with_rel_abundance <- calc_relative_abundance(qfeatures_obj, "protein")
#' 
#' # Access the relative abundance assay:
#' # rel_abundance <- assay(qf_with_rel_abundance, "protein_rel_abundance")
#' 
#' # Use with plotting functions:
#' # plot_relative_abundance(qf_with_rel_abundance, "protein_rel_abundance")
#' 
#' # Use with downstream analysis:
#' # perform_limma_analysis(qf_with_rel_abundance, "protein_rel_abundance", ~group, "groupB - groupA")
#'
#' @note
#' This function:
#' \itemize{
#'   \item Preserves NA values from the input data
#'   \item Normalizes each sample independently
#'   \item Handles zero values appropriately
#'   \item Maintains the original assay structure
#' }
#' For large datasets, consider:
#' \itemize{
#'   \item Filtering low-abundance features before calculation
#'   \item Using add_relative_abundance_assays() for multiple assays
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[add_relative_abundance_assays]{add_relative_abundance_assays}} for
#'     calculating relative abundances for multiple assays
#'   \item \code{\link[plot_relative_abundance]{plot_relative_abundance}} for visualizing
#'     the relative abundance data
#' }
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
