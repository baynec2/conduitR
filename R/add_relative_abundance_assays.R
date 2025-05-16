#' Add Relative Abundance Assays to QFeatures Object
#'
#' Adds new assays containing relative abundance values (as percentages) to a QFeatures object
#' for each specified taxonomic level. This function calculates the relative abundance of
#' features within each sample by normalizing the raw abundance values to percentages,
#' making it easier to compare the composition of features across samples.
#'
#' @param qf A QFeatures object containing the abundance data. The assays should contain
#'   non-negative abundance values (e.g., protein intensities or species counts).
#' @param assay_names A character vector specifying which assays to convert to relative
#'   abundances. By default, includes all taxonomic levels: domain, kingdom, phylum,
#'   class, order, family, genus, and species. Each assay should exist in the QFeatures
#'   object.
#'
#' @return A QFeatures object with additional assays named "{assay_name}_rel_abundance"
#'   for each input assay. Each new assay contains:
#'   \itemize{
#'     \item Values normalized to percentages (0-100%)
#'     \item Same dimensions as the input assay
#'     \item NA values preserved from the input
#'     \item Column sums equal to 100% for each sample (excluding NA values)
#'   }
#'
#' @export
#'
#' @examples
#' # Add relative abundance assays for all default taxonomic levels:
#' # qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj)
#' 
#' # Add relative abundance for specific assays:
#' # qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj,
#' #   assay_names = c("genus", "species"))
#' 
#' # The resulting assays can be used with plotting functions:
#' # plot_relative_abundance(qf_with_rel_abundance, "genus_rel_abundance")
#' 
#' # Or for downstream analysis:
#' # perform_limma_analysis(qf_with_rel_abundance, "species_rel_abundance", ~group, "groupB - groupA")
#'
#' @note
#' The input assays should contain non-negative abundance values. The function:
#' \itemize{
#'   \item Preserves NA values from the input data
#'   \item Handles each assay independently
#'   \item Normalizes each sample separately
#' }
#' For large datasets, consider processing one assay at a time to manage memory usage.
#'
#' @seealso \code{\link[plot_relative_abundance]{plot_relative_abundance}} for visualizing
#'   the relative abundance data
add_relative_abundance_assays <- function(qf, assay_names = c("domain",
                                                              "kingdom",
                                                              "phylum",
                                                              "class",
                                                              "order",
                                                              "family",
                                                              "genus",
                                                              "species")
                                          ){

  # Loop through all specified assay names, add relative abundance.
  for(i in assay_names){
    # Check if the specified assay exists in the QFeatures object
    if (!(i %in% names(qf))) {
      stop("The specified assay name does not exist in the QFeatures object.")
    }
    # Extract the assay data
    assay_data <- SummarizedExperiment::assay(qf, i)

    # Calculate column sums (total abundance for each sample)
    col_sums <- colSums(assay_data, na.rm = TRUE)

    # Normalize each value by its corresponding column sum to get relative abundance, multiply by 100 to get percentage
    rel_abundance <- QFeatures::sweep(assay_data, 2, col_sums, FUN = "/") * 100

    # Ensure the relative abundance is wrapped in a SummarizedExperiment
    rel_abundance_se <- SummarizedExperiment::SummarizedExperiment(assays = rel_abundance)

    # Add the relative abundance as a new assay in the QFeatures object
    qf <- QFeatures::addAssay(qf, rel_abundance_se, name = paste0(i, "_rel_abundance"))

  }

  # Return the modified QFeatures object
  return(qf)
}
