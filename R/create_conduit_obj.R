#' create_conduit_obj
#' create a conduit object.
#'
#' @param qfeatures.rds this is a qfeatures object (intended to be the output of
#' the conduit snakemake workflow)
#' @param metrics.tsv this is the metrics.tsv file created by the conduit snakemake workflow.
#' contains information of the number of proteins/ taxa (at all levels) that
#' were detected and
#' @param database_taxonomy.tsv this is the database taxonomy that was included in the
#' study.
#' @param detected_protein_metrics.tsv this is the proteins that were detected
#' @param precomputed_plots this is a list of precomputed plots (precomputed so
#' they don't take forever to render with conduit-GUI). Names of plots below.
#' * database_taxa_tree
#' * detected_taxa_tree
#' * difference_taxa_tree. Percentage relative to database.
#' * database_sunburst
#' * detected_sunburst
#' @returns
#' @export
#'
#' @examples
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
