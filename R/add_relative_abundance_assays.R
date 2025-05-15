#' add_relative_abundance_assays
#' add an assay containing relative abundance (as a proportion) to a qfeatures
#' object.
#'
#' @param qf
#' @param assay_name
#'
#' @returns
#' @export
#'
#' @examples
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
