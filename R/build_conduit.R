#' Build a Conduit Object from Multiple Data Sources
#'
#' Assembles a conduit object by combining a QFeatures object with protein taxonomy
#' and metrics data. This function integrates multiple data sources into a single,
#' organized data structure that maintains relationships between proteins, their
#' taxonomy, and associated metrics.
#'
#' @param QFeatures A QFeatures object containing proteomics data with protein
#'   abundances across samples.
#' @param diann_stats A named list of tibbles containing protein metrics (e.g.,
#'   identification scores, sequence coverage, etc.).
#' @param database A tibble containing protein IDs and the corresponding taxonomy ID.
#' @param annotations A tibble containing protein group, protein ID, annotations
#'   and taxonomy at the LCA level.
#' @param taxonomy A tibble containing the taxonomy detected in the experiment.
#' @param provenance Optional named list of provenance metadata created by
#'   [create_provenance()]. Defaults to NULL.
#'
#' @return A conduit object containing:
#'   \itemize{
#'     \item QFeatures: The proteomics data object
#'     \item metrics: Taxonomy information for all database proteins
#'     \item database: Protein metrics and quality scores
#'     \item annotations: annotation information for detected proteins
#'     \item taxonomy: taxonomy detected in the experiment
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' conduit_obj <- build_conduit_obj(
#'   QFeatures = qf_rds,
#'   diann_stats = "diann_stats.tsv",
#'   database = "database.tsv",
#'   annotations = "annotations.tsv",
#'   taxonomy = "taxonomy.tsv"
#' )
#' }
#'
#' @note
#' This function requires:
#' \itemize{
#'   \item A valid QFeatures object saved as .rds
#'   \item Tab-delimited text files for metrics and taxonomy
#'   \item Consistent protein identifiers across all input files
#' }
#' The function will read all files into memory, so ensure sufficient RAM is
#' available for large datasets.
build_conduit_obj <- function(QFeatures,
                              diann_stats,
                              database,
                              annotations,
                              taxonomy,
                              provenance = NULL) {

  # Create the  conduit object
  conduit_obj <- new("conduit",
    QFeatures = QFeatures,
    metrics = diann_stats,
    database = database,
    annotations = annotations,
    taxonomy = taxonomy,
    provenance = provenance
  )

  invisible(conduit_obj)
}
