#'Add Protein Coverage Metrics at All Taxonomic Levels
#'
#' Calculates the proportion of proteins detected for each taxonomic level
#' (domain, kingdom, phylum, class, order, family, genus, species) and adds
#' the results as separate metrics to the `metrics` list of the `conduit` object.
#'
#' @param conduit A `conduit` object containing proteomics data and database information.
#'
#' @return The original `conduit` object with new metrics added to
#'   \code{conduit@metrics}. Each metric is a tibble containing:
#'   \itemize{
#'     \item \code{taxon}: The taxon name at that level (e.g., species, genus)
#'     \item \code{n_proteins_db}: Total number of proteins in the database for that taxon
#'     \item \code{n_proteins_detected}: Number of detected proteins assigned to that taxon
#'     \item \code{coverage}: Percentage of proteins detected (\code{n_proteins_detected / n_proteins_db * 100})
#'   }
#'
#'   Metric names in the `metrics` list follow the pattern:
#'   \code{protein_coverage_<taxonomic_level>}, e.g., \code{protein_coverage_species}.
#'
#' @examples
#' \dontrun{
#' conduit_obj <- add_protein_coverage_all_taxa(conduit_obj)
#' head(conduit_obj@metrics$protein_coverage_species)
#' head(conduit_obj@metrics$protein_coverage_genus)
#' }
#'
#' @export
add_protein_coverage_taxa_metrics <- function(conduit) {

  taxonomic_levels <- c("domain", "kingdom", "phylum", "class",
                        "order", "family", "genus", "species")

  detected <- SummarizedExperiment::rowData(conduit@QFeatures[["protein_groups"]])[,c("Protein.Group",taxonomic_levels)] |>
    tibble::as_tibble()

  database <- conduit@database

  for (level in taxonomic_levels) {
    if (level %in% colnames(detected)) {

      # Total proteins per taxon in database
      db_count <- database |>
        dplyr::group_by(taxon = .data[[level]]) |>
        dplyr::summarise(n_proteins_db = dplyr::n_distinct(protein_id),
                         .groups = "drop")

      # Detected proteins per taxon
      detected_count <- detected |>
        tidyr::separate_rows(Protein.Group, sep = ";") |>
        dplyr::group_by(taxon = .data[[level]]) |>
        dplyr::summarise(n_proteins_detected = dplyr::n_distinct(Protein.Group),
                         .groups = "drop")

      # Combine and calculate coverage
      coverage_tbl <- dplyr::left_join(db_count, detected_count, by = "taxon") |>
        dplyr::mutate(coverage = round(n_proteins_detected / n_proteins_db * 100, 2))

      # Add to metrics list
      metric_name <- paste0("protein_coverage_", level)
      conduit@metrics[[metric_name]] <- coverage_tbl
    }
  }

  return(conduit)
}
