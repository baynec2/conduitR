#' Count Proteins per Taxonomic Group
#'
#' Counts the number of proteins present in each taxonomic group across all
#' taxonomic levels. This function is used internally by other functions to
#' calculate protein detection rates and create taxonomic visualizations.
#'
#' @param tibble A data frame or tibble containing protein taxonomy information.
#'   Must have the following columns:
#'   \itemize{
#'     \item organism_type: Type of organism (e.g., "Bacteria", "Archaea")
#'     \item domain: Taxonomic domain
#'     \item kingdom: Taxonomic kingdom
#'     \item phylum: Taxonomic phylum
#'     \item class: Taxonomic class
#'     \item order: Taxonomic order
#'     \item family: Taxonomic family
#'     \item genus: Taxonomic genus
#'     \item species: Taxonomic species
#'   }
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item All taxonomic columns from the input
#'     \item n: Number of proteins in each taxonomic group
#'   }
#'   The result is grouped by all taxonomic levels, allowing for hierarchical
#'   analysis of protein counts.
#'
#' @export
#'
#' @examples
#' # Count proteins in each taxonomic group:
#' # counts <- count_proteins_per_taxon(taxonomy_data)
#' 
#' # The results can be used with other functions:
#' # plot_percent_detected_taxa_tree(conduit_obj)
#' # plot_sunburst(taxonomy_data)
count_proteins_per_taxon = function(tibble){

  sum = tibble |>
    dplyr::group_by(.data$organism_type,.data$domain,.data$kingdom,
                    .data$phylum,.data$class,.data$order,.data$family,
                    .data$genus,.data$species) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::ungroup()

  return(sum)

}
