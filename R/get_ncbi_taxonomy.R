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
#'     \item species: Scientific name at *species* rank
#'     \item strain: Infraspecific designation, or NA when the queried taxon is
#'       itself at species rank
#'     \item rank: Taxonomic rank (e.g., domain, kingdom, phylum)
#'     \item name: Scientific name at each taxonomic level
#'   }
#'   The data frame includes all taxonomic levels from domain to genus,
#'   with missing ranks filled as NA.
#'
#' @details
#' `species` always holds the name at species rank, never a strain designation.
#' Many UniProt reference proteomes are registered under strain taxids (e.g.
#' 349741, `Akkermansia muciniphila ATCC BAA-835`, rank `strain`), whose own
#' scientific name carries the strain. Writing that name into `species` splits
#' one species into several, so a peptide shared between two strains has no
#' common species and its LCA falls to genus — dropping it from every
#' species-level aggregation.
#'
#' The species-rank ancestor is therefore resolved from the record's own `Rank`
#' plus its `LineageEx`, both already present in the fetched XML:
#' \itemize{
#'   \item queried taxon is rank `species` → `species` is its own name,
#'     `strain` is NA;
#'   \item queried taxon is below species and its lineage carries a species-rank
#'     ancestor → `species` is that ancestor, `strain` is the infraspecific
#'     remainder;
#'   \item otherwise → the record's own name is kept in `species` unchanged,
#'     since no better information is available.
#' }
#'
#' Because the rule keys on rank rather than on names, informal designations are
#' unaffected: `Lachnospiraceae bacterium A2` (taxid 397290) is itself rank
#' `species` with no species-rank ancestor, so it is never folded into
#' `Lachnospiraceae bacterium` (taxid 1898203). Both are left exactly as NCBI
#' reports them.
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

        # Extract scientific name of the queried organism. NOTE: this is the
        # record's own name, which for a strain taxid carries the strain
        # designation -- it is not necessarily a species-rank name.
        own_name <- XML::xpathSApply(
          xml_parsed,
          "//TaxaSet/Taxon/ScientificName",
          XML::xmlValue
        )[1]

        # Rank of the queried organism itself (e.g. "species", "strain").
        own_rank <- XML::xpathSApply(
          xml_parsed,
          "//TaxaSet/Taxon/Rank",
          XML::xmlValue
        )[1]

        # Extract tax ID
        organism_id <- XML::xpathSApply(
          xml_parsed,
          "//TaxaSet/Taxon/TaxId",
          XML::xmlValue
        )[1]

        # Resolve the species-rank name and any infraspecific remainder. See
        # @details: keyed on rank, so informal names are never merged.
        lineage_species <- names[ranks == "species"]
        lineage_species <- if (length(lineage_species) > 0) {
          lineage_species[1]
        } else {
          NA_character_
        }

        if (identical(own_rank, "species")) {
          species <- own_name
          strain <- NA_character_
        } else if (!is.na(lineage_species)) {
          species <- lineage_species
          # Strip the species name to leave just the infraspecific part; keep
          # the full name when it is not a clean prefix (NCBI is not uniform).
          prefix <- paste0(lineage_species, " ")
          strain <- if (startsWith(own_name, prefix)) {
            trimws(substring(own_name, nchar(prefix) + 1L))
          } else {
            own_name
          }
        } else {
          # Below species but no species-rank ancestor available: preserve the
          # record's own name rather than discarding it.
          species <- own_name
          strain <- NA_character_
        }

        # Combine all results into one dataframe
        # Define expected taxonomy ranks
        expected_ranks <- c(
          "domain", "kingdom", "phylum", "class", "order",
          "family", "genus"
        )

        taxonomy_df <- data.frame(
          organism_id = organism_id,
          species = species,
          strain = strain,
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
            strain = strain,
            name = NA_character_
          ))

        # Index by character: a numeric taxid would extend the list to that
        # length, allocating millions of empty slots.
        taxonomy_list[[as.character(id)]] <- taxonomy_df
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
      "order", "family", "genus", "species", "strain"
    ))

  return(final_df)
}
