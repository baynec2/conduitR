#' Plot Taxonomic Tree with Detection Rates
#'
#' Creates an interactive heat tree visualization showing the taxonomic distribution
#' of detected proteins, with node sizes and colors representing the percentage of
#' proteins detected at each taxonomic level. This helps visualize the coverage
#' of protein detection across different taxonomic groups.
#'
#' @param conduit_obj A conduit object containing protein taxonomy information
#' @param type Character string specifying how to calculate detection rates:
#'   \itemize{
#'     \item "multiple_proteins_in_group": Count proteins that appear in multiple
#'       samples within a group (default)
#'     \item "any_protein_in_group": Count any protein detected in a group
#'   }
#' @param ... Additional arguments passed to metacoder::heat_tree()
#'
#' @return An interactive heat tree plot where:
#'   \itemize{
#'     \item Nodes represent taxonomic groups
#'     \item Node size indicates the percentage of proteins detected
#'     \item Node color indicates the percentage of proteins detected
#'     \item Node labels show taxonomic names
#'     \item The tree structure shows taxonomic relationships
#'   }
#'
#' @export
#'
#' @examples
#' # Basic taxonomic tree with default settings:
#' # plot_percent_detected_taxa_tree(conduit_obj)
#' 
#' # Using alternative detection rate calculation:
#' # plot_percent_detected_taxa_tree(conduit_obj, type = "any_protein_in_group")
#' 
#' # Customizing the heat tree appearance:
#' # plot_percent_detected_taxa_tree(conduit_obj,
#' #                                node_label = "taxon_names",
#' #                                node_size_range = c(0.01, 0.1),
#' #                                node_color_range = c("white", "red"))
plot_percent_detected_taxa_tree = function(conduit_obj,
                                           type ="multiple_proteins_in_group",
                                           ...){

  # Accessing data stored in the conduit object
  database_protein_taxonomy = slot(conduit_obj,"database_protein_taxonomy")
  detected_protein_taxonomy = slot(conduit_obj,"detected_protein_taxonomy")

  # Calculating the percent of proteins detected.
  pd = calc_percent_proteins_detected(conduit_obj,type = type)

  # Create taxonomy string for use with metacoder
  taxonomy_data <- pd |>
    dplyr::mutate(taxonomy = paste0(
      "organism_type__", organism_type, ";",
      "domain__", domain, ";",
      "kingdom__", kingdom, ";",
      "phylum__", phylum, ";",
      "class__", class, ";",
      "order__", order, ";",
      "family__", family, ";",
      "genus__", genus, ";",
      "species__", species
    )) |>
    dplyr::select(species, taxonomy, n_detected,n_in_db,percent_detected)

  # Parse the taxonomy string to create the taxmap object
  taxmap <- metacoder::parse_tax_data(
    taxonomy_data,
    class_cols = "taxonomy",  # Column containing the concatenated taxonomy string
    class_sep = ";",  # Separator used in taxonomy string
    class_key = c(
      tax_rank = "taxon_rank",  # Key to capture the rank
      tax_name = "taxon_name"   # Key to capture the taxon name
    ),
    class_regex = "^(.+)__(.+)$"  # Regex to separate rank and name
  )

  # Extracting the amount in the db
  taxmap$data$tax_abund_db <- metacoder::calc_taxon_abund(taxmap, "tax_data",
                                               cols = "n_in_db",
                                               groups = "n_in_db")
  # Getting the amount detected
  taxmap$data$tax_abund_total <- metacoder::calc_taxon_abund(taxmap, "tax_data",
                                                  cols = "n_detected",
                                                  groups = "total_count")
  # Calculating the percent detected
 taxmap$data$tax_abund_per <- tibble::tibble(taxon_id = taxmap$data$tax_abund_total$taxon_id,
                                             percent_proteins_detected = (taxmap$data$tax_abund_total$total_count / taxmap$data$tax_abund_db$n_in_db) * 100)

  # plot the heat tree
 p1 = metacoder::heat_tree(taxmap,
                node_label = taxon_names,
                node_size = percent_proteins_detected,
                node_color = percent_proteins_detected,
                ...)

 return(p1)
}

