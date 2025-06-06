#' Create a Conduit Object from Multiple Data Sources
#'
#' Assembles a conduit object by combining a QFeatures object with protein taxonomy
#' and metrics data. This function integrates multiple data sources into a single,
#' organized data structure that maintains relationships between proteins, their
#' taxonomy, and associated metrics.
#'
#' @param QFeatures.rds Character string specifying the file path to a saved QFeatures
#'   object (.rds file). The QFeatures object should contain proteomics data with
#'   protein abundances across samples.
#' @param combined_metrics.tsv Character string specifying the file path to a
#'   tab-delimited text file containing protein metrics (e.g., identification scores,
#'   sequence coverage, etc.).
#' @param database_protein_taxonomy.tsv Character string specifying the file path to
#'   a tab-delimited text file containing taxonomy information for all proteins in
#'   the database.
#' @param detected_protein_taxonomy.tsv Character string specifying the file path to
#'   a tab-delimited text file containing taxonomy information for detected proteins.
#'
#' @return A conduit object containing:
#'   \itemize{
#'     \item QFeatures: The proteomics data object
#'     \item database_protein_taxonomy: Taxonomy information for all database proteins
#'     \item combined_metrics: Protein metrics and quality scores
#'     \item detected_protein_taxonomy: Taxonomy information for detected proteins
#'   }
#'
#' @export
#'
#' @examples
#' # Create a conduit object from saved files:
#' # conduit_obj <- create_conduit_obj(
#' #   QFeatures.rds = "proteomics_data.rds",
#' #   combined_metrics.tsv = "protein_metrics.tsv",
#' #   database_protein_taxonomy.tsv = "all_proteins_taxonomy.tsv",
#' #   detected_protein_taxonomy.tsv = "detected_proteins_taxonomy.tsv"
#' # )
#'
#' # The resulting object can be used with other conduitR functions:
#' # plot_taxa_tree(conduit_obj)
#' # plot_relative_abundance(conduit_obj, "protein")
#'
#' @note
#' This function requires:
#' \itemize{
#'   \item A valid QFeatures object saved as .rds
#'   \item Tab-delimited text files for metrics and taxonomy
#'   \item Consistent protein identifiers across all input files
#' }
#' The function will read all files into memory, so ensure sufficient RAM is available
#' for large datasets.
create_conduit_obj <- function(QFeatures.rds,
                               combined_metrics.tsv,
                               database_protein_taxonomy.tsv,
                               detected_protein_taxonomy.tsv) {
  # Reading in the files
  QFeatures <- readRDS(QFeatures.rds)
  database_protein_taxonomy <- readr::read_tsv(database_protein_taxonomy.tsv)
  combined_metrics <- readr::read_tsv(combined_metrics.tsv)
  detected_protein_taxonomy <- readr::read_tsv(detected_protein_taxonomy.tsv)

  # Create the  conduit object
  conduit_obj <- new("conduit",
    QFeatures = QFeatures,
    database_protein_taxonomy = database_protein_taxonomy,
    combined_metrics = combined_metrics,
    detected_protein_taxonomy = detected_protein_taxonomy
  )
  return(conduit_obj)
}
