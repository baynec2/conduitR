#' Calculate Protein Detection Rates by Taxonomy
#'
#' Computes the percentage of proteins detected at each taxonomic level by comparing
#' the number of detected proteins against the total number of proteins in the database
#' for each taxonomic group.
#'
#' @param conduit_obj A conduit object containing:
#'   \itemize{
#'     \item database_protein_taxonomy: Taxonomy information for all proteins in the database
#'     \item detected_protein_taxonomy: Taxonomy information for detected proteins
#'     \item QFeatures: Protein group information (for "one_protein_in_group" type)
#'   }
#' @param type Character string specifying how to count detected proteins:
#'   \itemize{
#'     \item "multiple_proteins_in_group": Count all proteins in detected protein groups (default)
#'     \item "one_protein_in_group": Count only uniquely identified proteins (no shared peptides)
#'   }
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item Taxonomic columns (organism_type, domain, kingdom, phylum, class, order, family, genus, species)
#'     \item n_in_db: Number of proteins in the database for each taxon
#'     \item n_detected: Number of proteins detected for each taxon
#'     \item percent_detected: Percentage of proteins detected (n_detected/n_in_db * 100)
#'   }
#'
#' @export
#'
#' @examples
#' # Calculate detection rates using all proteins in groups:
#' # rates <- calc_percent_proteins_detected(conduit_obj)
#' 
#' # Calculate detection rates using only uniquely identified proteins:
#' # rates <- calc_percent_proteins_detected(conduit_obj, type = "one_protein_in_group")
#' 
#' # The results can be used with plot_percent_detected_taxa_tree():
#' # plot_percent_detected_taxa_tree(conduit_obj)
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

