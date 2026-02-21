#' Fetch Complete Taxonomic Information from NCBI
#'
#' Retrieves comprehensive taxonomic information for a list of NCBI taxonomy IDs
#' using the NCBI Entrez API. This function is useful for obtaining complete
#' taxonomic lineages and scientific names for organisms.
#'
#' @param ncbi_ids A numeric vector of NCBI taxonomy IDs to query
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item organism_id: NCBI taxonomy ID
#'     \item species: Scientific name of the organism
#'     \item rank: Taxonomic rank (e.g., domain, kingdom, phylum)
#'     \item name: Scientific name at each taxonomic level
#'   }
#'   The data frame includes all taxonomic levels from domain to genus,
#'   with missing ranks filled as NA.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single organism (human)
#' taxonomy <- get_ncbi_taxonomy(9606)
#'
#' # Multiple organisms (human and mouse)
#' taxonomy <- get_ncbi_taxonomy(c(9606, 10090))
#'
#' # Use with plotting functions
#' plot_taxa_tree(taxonomy)
#' plot_sunburst(taxonomy)
#' }
#' # Requires internet; NCBI may rate-limit requests.
#'
get_ncbi_taxonomy <- function(ncbi_ids) {
  taxonomy_list <- list()
  for (id in ncbi_ids) {
    tryCatch(
      {
        # Query NCBI Taxonomy
        xml_data <- rentrez::entrez_fetch(db = "taxonomy",
                                          id = id,
                                          rettype = "xml")

        # Parse the XML data using XML package
        xml_parsed <- XML::xmlParse(xml_data)

        # Extract taxonomic lineage
        ranks <- XML::xpathSApply(
          xml_parsed,
          "//LineageEx/Taxon/Rank",
          XML::xmlValue
        )
        names <- XML::xpathSApply(
          xml_parsed,
          "//LineageEx/Taxon/ScientificName",
          XML::xmlValue
        )

        # Extract scientific name of the queried organism
        species <- XML::xpathSApply(
          xml_parsed,
          "//TaxaSet/Taxon/ScientificName",
          XML::xmlValue
        )

        # Extract tax ID
        organism_id <- XML::xpathSApply(
          xml_parsed,
          "//TaxaSet/Taxon/TaxId",
          XML::xmlValue
        )

        # Combine all results into one dataframe
        # Define expected taxonomy ranks
        expected_ranks <- c(
          "domain", "kingdom", "phylum", "class", "order",
          "family", "genus"
        )

        taxonomy_df <- data.frame(
          organism_id = organism_id,
          species = species,
          rank = ranks,
          name = names,
          stringsAsFactors = FALSE
        ) |>
          # Ensure ranks are only those in expected_ranks
          dplyr::mutate(rank = ifelse(rank %in% expected_ranks, rank, NA)) |>
          tidyr::drop_na(rank) |> # Remove any entries with NA rank
          # Ensure all expected ranks are present, filling missing ones with NA
          tidyr::complete(rank = expected_ranks, fill = list(
            organism_id = organism_id,
            species = species,
            name = NA_character_
          ))

        taxonomy_list[[id]] <- taxonomy_df
      },
      error = function(e) {
        log_with_timestamp(paste(
          "Error fetching taxonomy for ID:", id, "|", conditionMessage(e)
        ))      }
    )
  }

  final_df <- dplyr::bind_rows(taxonomy_list) |>
    tidyr::pivot_wider(names_from = rank, values_from = name) |>
    dplyr::select(c(
      "organism_id", "domain", "kingdom", "phylum", "class",
      "order", "family", "genus", "species"
    ))

  return(final_df)
}
