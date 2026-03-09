#' Add Relative Abundance Assays for Multiple Taxonomic Levels
#'
#' Wrapper around \code{\link{add_relative_abundance_assay}} that converts all
#' standard taxonomic-level assays present in a QFeatures object to relative
#' abundance (percentage) assays.
#'
#' @param qf A QFeatures object.
#' @param assay_names Character vector of assay names to convert. Defaults to all
#'   standard taxonomic levels that are present in \code{qf}.
#' @return A QFeatures object with additional \code{{assay_name}_rel_abundance} assays.
#' @export
#'
#' @examples
#' # qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj)
#' # qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj,
#' #   assay_names = c("genus", "species"))
add_relative_abundance_assays <- function(qf,
                                          assay_names = c("domain", "kingdom", "phylum",
                                                          "class", "order", "family",
                                                          "genus", "species")) {
  assay_names <- intersect(assay_names, names(qf))
  for (assay_name in assay_names) {
    qf <- add_relative_abundance_assay(qf, assay_name = assay_name)
  }
  return(qf)
}

#' Add Relative Abundance Assay for a Single Taxonomic Level
#'
#' Adds a new assay containing relative abundance values (as percentages) for
#' one taxonomic level in a QFeatures object.
#'
#' @param qf A QFeatures object.
#' @param assay_name Name of the taxonomic assay to convert. One of: "domain",
#'   "kingdom", "phylum", "class", "order", "family", "genus", "species".
#' @return A QFeatures object with an additional \code{{assay_name}_rel_abundance} assay.
#' @export
add_relative_abundance_assay <- function(qf, assay_name = "phylum") {
  allowed_assays <- c(
    "domain", "kingdom", "phylum", "class", "order", "family",
    "genus", "species"
  )

  if (!(assay_name %in% allowed_assays)) {
    stop("The specified assay name does not exist in the QFeatures object.")
  }

  # Extract the assay
  se <- qf[[assay_name]]
  assay_data <- SummarizedExperiment::assay(se)

  # Compute relative abundance (%)
  col_sums <- colSums(assay_data, na.rm = TRUE)
  rel_abundance <- sweep(assay_data, 2, col_sums, FUN = "/") * 100

  # Create a new SummarizedExperiment with original metadata
  rel_abundance_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(rel_abundance = rel_abundance),
    rowData = SummarizedExperiment::rowData(se),
    colData = SummarizedExperiment::colData(se)
  )


  # Add the relative abundance as a new assay in the QFeatures object
  qf <- QFeatures::addAssay(qf, rel_abundance_se, name = paste0(
    assay_name,
    "_rel_abundance"
  ))

  # Return the modified QFeatures object
  return(qf)
}
