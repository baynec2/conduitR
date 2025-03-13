#' create_probiomeconduit_obj
#' create a pbiomeconduit obj
#'
#' @param qfeature.rda
#' @param database_taxonomy.tsv
#' @param database_metrics.tsv
#' @param detected_protein_taxonomy.tsv
#' @param detected_protein_metrics.tsv
#'
#' @returns
#' @export
#'
#' @examples
create_probiomeconduit_obj = function(qfeature.rds,
                                   database_taxonomy.tsv,
                                   database_metrics.tsv,
                                   detected_protein_taxonomy.tsv,
                                   detected_protein_metrics.tsv){
  #Reading in the files
  QFeatures = readRDS(qfeature.rds)
  database_taxonomy = readr::read_tsv(database_taxonomy.tsv)
  database_metrics = readr::read_tsv(database_metrics.tsv)
  detected_protein_taxonomy = readr::read_tsv(detected_protein_taxonomy.tsv)
  detected_protein_metrics = readr::read_tsv(detected_protein_metrics.tsv)

# Create the probiotics conduit object
probiome_conduit_obj <- new("probiomeconduit",
                            QFeatures = QFeatures,
                            database_taxonomy = database_taxonomy,
                            database_metrics = database_metrics,
                            detected_protein_taxonomy = detected_protein_taxonomy,
                            detected_protein_metrics = detected_protein_metrics)

return(probiome_conduit_obj)
}
