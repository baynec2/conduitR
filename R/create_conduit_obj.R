#' Create a Conduit Object
#'
#' Creates a conduit object that combines proteomics data, taxonomy information,
#' and analysis metrics into a single object for downstream analysis and visualization.
#' This object serves as the main data structure for the conduit package.
#'
#' @param QFeatures.rds Character string specifying the path to an RDS file containing
#'   a QFeatures object with proteomics data (typically output from the conduit
#'   snakemake workflow)
#' @param combined_metrics.tsv Character string specifying the path to a TSV file
#'   containing metrics about detected proteins and taxa at all taxonomic levels
#' @param database_protein_taxonomy.tsv Character string specifying the path to a TSV file
#'   containing taxonomy information for all proteins in the database
#' @param detected_protein_taxonomy.tsv Character string specifying the path to a TSV file
#'   containing taxonomy information for detected proteins
#'
#' @return A conduit object containing:
#'   \itemize{
#'     \item QFeatures: A QFeatures object with proteomics data
#'     \item database_protein_taxonomy: A data frame with taxonomy information for all database proteins
#'     \item combined_metrics: A data frame with detection metrics
#'     \item detected_protein_taxonomy: A data frame with taxonomy information for detected proteins
#'   }
#'
#' @export
#'
#' @examples
#' # Create a conduit object from workflow outputs:
#' # conduit_obj <- create_conduit_obj(
#' #   QFeatures.rds = "results/qfeatures.rds",
#' #   combined_metrics.tsv = "results/metrics.tsv",
#' #   database_protein_taxonomy.tsv = "results/database_taxonomy.tsv",
#' #   detected_protein_taxonomy.tsv = "results/detected_taxonomy.tsv"
#' # )
#' 
#' # The resulting object can be used with various analysis and plotting functions:
#' # - plot_percent_detected_taxa_tree()
#' # - plot_taxa_tree()
#' # - plot_sunburst()
#' # - calc_percent_proteins_detected()
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
