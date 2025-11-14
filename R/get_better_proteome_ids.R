#' Select the best UniProt proteome for a given taxonomic identifier.
#'
#' Sometimes, the proteome IDs returned for a specific NCBI taxonomic identifier
#' are not available in UniProtKB. To maximize the information obtained from
#' these identifiers, this function applies the following selection logic:
#'
#' 1. **Strain-level taxonomic IDs:**
#'    - If a proteome exists that is classified as "Representative proteome",
#'      "Reference and representative proteome", or "Other proteome", it will be
#'      used.
#'    - If no such proteome exists, the function moves up to the species level
#'      and selects proteomes in the following order of preference:
#'      "Representative proteome", "Reference and representative proteome",
#'      "Other proteome", "Redundant proteome", and "Excluded proteome".
#'
#' 2. **Species-level taxonomic IDs:**
#'    - If a redundant proteome exists, it will be used as the selected
#'      proteome.
#'
#' The first proteome meeting these criteria will be returned as `selected_proteome_id`.
#'
#' @param proteome_id_df Data frame containing annotations produced by the
#'        `get_proteome_ids()` function.
#' @param parallel Logical; whether to perform selection in parallel (default: FALSE).
#'
#' @return A data frame with an additional column `selected_proteome_id` indicating
#'         the chosen proteome for each taxonomic identifier.
#'
#' @export
#'
#' @examples
#' # Data format compatible with get_better_proteome_ids()
#' proteome_id_df <- tibble::tibble(
#'   proteome_id = "UP000027861",
#'   organism_id = 1339341,
#'   organism = "Parabacteroides distasonis str. 3776 Po2 i",
#'   protein_count = 4580,
#'   proteome_type = "Redundant proteome",
#'   redundant_to = "UP000027850",
#'   genome_assembly_id = "GCA_000699745.1",
#'   genome_assembly_level = "full",
#'   annotation_score = NA)
#'
#' # Example function call
#' better_proteome_ids <- get_better_proteome_ids(proteome_id_df)

get_better_proteome_ids <- function(proteome_id_df,
                                    parallel = FALSE) {

  # Extract organism ids
  organism_ids <- proteome_id_df$organism_id

  log_with_timestamp("Getting parent taxonomy for each proteome")

  # Get taxonomy info for each organism
  taxon_df <- purrr::map_dfr(organism_ids, get_parent_taxonomy_id, .progress = TRUE)

  # Merge proteome info with taxonomy
  proteome_taxa <- dplyr::left_join(
    proteome_id_df,
    taxon_df,
    by = c("organism_id" = "child_id")
  ) |>
    dplyr::mutate(selected_proteome_id = dplyr::coalesce(redundant_to,
                                                         proteome_id))

  # ----------------------------
  # Strain-level handling.
  # ----------------------------
  strain_df <- proteome_taxa |> dplyr::filter(child_rank == "strain")

  if (nrow(strain_df) > 0) {

    good_pt_filter <- c("Representative proteome",
                        "Reference and representative proteome",
                        "Other proteome")

    log_with_timestamp("Determining what strains have adequate proteomes")

    # Good strains
    good_strains <- strain_df |> dplyr::filter(proteome_type %in% good_pt_filter)

    # Bad strains (non-representative)
    bad_strains <- strain_df |> dplyr::filter(proteome_type %!in% good_pt_filter)

    # Get parent species proteomes for bad strains
    parents_of_bad <- unique(bad_strains$parent_id)

    log_with_timestamp("Getting species level proteomes for strains with bad uniprot proteomes")


    if (length(parents_of_bad) > 0) {
      parent_proteomes <- get_proteome_ids_from_organism_ids(parents_of_bad, parallel = parallel) |>
        dplyr::mutate(
          selected_proteome_id = dplyr::coalesce(redundant_to, proteome_id),
          child_rank = "species"
        )

      good_parents <- parent_proteomes |> dplyr::filter(proteome_type %in% good_pt_filter)
    } else {
      parent_proteomes <- dplyr::tibble()
      good_parents <- dplyr::tibble()
    }

    # Strains whose parent species do NOT have good proteomes
    strain_level_okay <- bad_strains |>
      dplyr::filter(parent_id %!in% good_parents$organism_id) |>
      dplyr::mutate(selected_proteome_id = dplyr::coalesce(redundant_to, proteome_id))

    # Combine all resolved strain-level IDs
    strain_resolved <- dplyr::bind_rows(
      good_strains,
      good_parents,
      strain_level_okay
    )
  } else {
    strain_resolved <- dplyr::tibble()
  }

  # ----------------------------
  # Non-strain level (species/subspecies)
  # ----------------------------
  species_df <- proteome_taxa |> dplyr::filter(child_rank != "strain") |>
    dplyr::filter(!is.na(selected_proteome_id))

  # Combine all IDs
  resolved_df <- dplyr::bind_rows(strain_resolved, species_df) |>
    dplyr::distinct(selected_proteome_id, .keep_all = TRUE)

  log_with_timestamp(sprintf("Returning %d selected proteome IDs", nrow(resolved_df)))

  return(resolved_df)
}
