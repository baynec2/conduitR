#' Normalize assay abundance to species abundance
#'
#' It may be adventageous to normalize protein intensities to total organism
#' abundance, paticularly when assessing a question of the type
# “Is the expression of geneA from SPECIESX higher in SPECIESX in 92 COMMUNITY1
#` as compared to SPECIESX in COMMUNITY2?” or “Which genes differ in expression
# `93 between COMMUNITY1 and COMMUNITY2on the species level?”,
#' . See https://peerj.com/preprints/2846/#.
#'
#' @param qf QFeatures
#' @param protein_assay assay holding log2 protein counts
#' @param species_assay assay holding log2 species count
#' @param species_col column in rowData of protein_assay indicating species
#' the feature was derived from.
#' @param out_assay name of assay normalized to species abundance.
#' @return The input `QFeatures` object with an additional assay (name given by
#'   `out_assay`) containing protein abundances normalized to species abundance.
#' @export
#'
#' @examples
#' \dontrun{
#' qf <- normalize_assay_to_species(qf,
#'   protein_assay = "protein_groups_log2",
#'   species_assay = "species_log2",
#'   out_assay = "protein_groups_log2_species_normalized"
#' )
#' }
normalize_assay_to_species <- function(qf,
                                          protein_assay = "protein_groups_log2",
                                          species_assay = "species_log2",
                                          species_col = "species",
                                          out_assay = "protein_groups_log2_species_normalized") {

  # Basic validation
  if(!protein_assay %in% names(qf) ){
    stop("indicated protein assay not present in supplied QFeatures")
  }
  if(!species_assay %in% names(qf)){
    stop("indicated species assay not present in supplied QFeatures")
  }

  # Extract matrices
  pg <- SummarizedExperiment::assay(qf, protein_assay)
  sp <- SummarizedExperiment::assay(qf, species_assay)

  # Species mapping
  pg_species <- SummarizedExperiment::rowData(qf[[protein_assay]])[[species_col]]
  if (is.null(pg_species)) stop("Species column not found in rowData")

  # Check that non-NA species exist in species assay
  missing_species <- setdiff(unique(pg_species[!is.na(pg_species)]),
                             rownames(sp))

  if (length(missing_species) > 0) stop(
    "These species are missing from species assay: ",
    paste(missing_species, collapse = ", ")
  )

  # Initialize normalized matrix with NA
  pg_norm <- matrix(
    NA_real_,
    nrow = nrow(pg),
    ncol = ncol(pg),
    dimnames = dimnames(pg)
  )

  # Only normalize rows with valid species
  valid <- !is.na(pg_species)
  sp_expanded <- sp[pg_species[valid], , drop = FALSE]
  pg_norm[valid, ] <- pg[valid, ] - sp_expanded

  # Wrap in a SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(norm = pg_norm),
    rowData = rowData(qf[[protein_assay]]),
    colData = colData(qf)
  )

  # Add back to QFeatures
  qf <- QFeatures::addAssay(qf, se, name = out_assay)

  return(qf)
}
