#' get_proteome_ids_from_organism_ids
#'
#' @param organism_ids
#' @param batch_size
#' @param parallel
#'
#' @returns
#' @export
#'
#' @examples
get_proteome_ids_from_organism_ids <- function(organism_ids,
                                               parallel = TRUE) {

  if (parallel) {
    future::plan(future::multisession, workers = future::availableCores() - 1)
    results <- furrr::future_map_dfr(organism_ids, get_proteome_id_from_organism_id,
      .progress = TRUE
    )
    future::plan(future::sequential) # Reset to sequential processing
  } else {
    results <- purrr::map_dfr(organism_ids, get_proteome_id_from_organism_id,
      .progress = TRUE
    )
  }

  # Getting list of all reference proteome, pull reference if it exists, if not
  # take largest
  cat("Getting Reference Proteomes")
  reference <- get_all_reference_proteomes() |> dplyr::pull(Proteome_ID)
  cat("Reference Proteomes Sucessfully Retrieved")

  # Annotate with whether it is a reference proteome or not.
  results <- results |>
    dplyr::mutate(reference = dplyr::case_when(
      `Proteome Id` %in% reference ~ TRUE,
      TRUE ~ FALSE
    )) |>
    dplyr::group_by(`Organism Id`) |>
    dplyr::filter(if (any(reference)) reference else `Protein count` == max(`Protein count`)) |>
    dplyr::ungroup()

  return(results)
}
