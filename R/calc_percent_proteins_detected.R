#' calc_percent_proteins_detected
#'
#' calculate the percent of proteins detected.
#'
#' @param conduit_obj tibble with columns named...
#' @param type multiple_proteins_in_group", "one_protein_in_group
#'
#' @returns
#' @export
#'
#' @examples
calc_percent_proteins_detected = function(conduit_obj,
                                          type = "multiple_proteins_in_group"){

  if (type == "multiple_proteins_in_group") {
    # Accessing data stored in the conduit object
    database_protein_taxonomy <- slot(conduit_obj, "database_protein_taxonomy")
    detected_protein_taxonomy <- slot(conduit_obj, "detected_protein_taxonomy")
  } else if (type == "one_protein_in_group") {
    # Finding unique protein groups
    pg <- row.names(SummarizedExperiment::assay(slot(conduit_obj, "QFeatures"), "protein_group"))
    uniquely_ided <- pg[!grepl(".*;.*", pg)]  # Unique ided pg don't have a ;

    # Accessing data stored in the conduit object
    database_protein_taxonomy <- slot(conduit_obj, "database_protein_taxonomy")

    detected_protein_taxonomy <- slot(conduit_obj, "detected_protein_taxonomy") |>
      dplyr::filter(protein_id %in% uniquely_ided)
  } else {
    stop("type must be either multiple_proteins_in_group or one_protein_in_group")
  }

  # Read in count of proteins per taxa that were in the database
  db = count_proteins_per_taxon(database_protein_taxonomy) |>
    dplyr::rename(n_in_db =n)

  # Read in count of proteins per taxa that were detected
  detected = count_proteins_per_taxon(detected_protein_taxonomy) |>
    dplyr::rename(n_detected = n)

  # Combining db and detected into one df. Calculating percent detected
  all = dplyr::left_join(db,detected, by = c("organism_type",
                                             "domain",
                                             "kingdom",
                                             "phylum",
                                             "class",
                                             "order",
                                             "family",
                                             "genus",
                                             "species")) |>
    dplyr::mutate(percent_detected = n_detected/n_in_db * 100)

  return(all)

}

