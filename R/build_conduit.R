#' Build a Conduit Object from Multiple Data Sources
#'
#' Assembles a conduit object by combining a QFeatures object with protein taxonomy
#' and metrics data. This function integrates multiple data sources into a single,
#' organized data structure that maintains relationships between proteins, their
#' taxonomy, and associated metrics.
#'
#' @param QFeatures  Character string specifying the file path to a saved QFeatures
#'   object (.rds file). The QFeatures object should contain proteomics data with
#'   protein abundances across samples.
#' @param diann_stats Character string specifying the file path to a
#'   tab-delimited text file containing protein metrics (e.g., identification scores,
#'   sequence coverage, etc.).
#' @param database_fp Character string specifying the file path to
#'   a tab-delimited text file containing protein ids and the corresponding taxonomy id.
#' @param annotations_fp Character string specifying the file path to
#'   a tab-delimited text file containing protein group, protein id, annnotations
#'   and taxonomy at the lca level/].
#' @param taxonomy_fp Character string specifying the file path to a
#' tab-delimited text file containing the taxonomy detected in the experiment
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
                              taxonomy) {

  # Create the  conduit object
  conduit_obj <- new("conduit",
    QFeatures = QFeatures,
    metrics = metrics,
    database = database,
    annotations = annotations,
    taxonomy = taxonomy
  )

  invisible(conduit_obj)
}
