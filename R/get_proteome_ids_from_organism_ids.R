#' get_proteome_ids_from_organism_ids
#'
#' get proteome ids from a vector of organism ids. This function will return a
#' reference proteome if one is available, else it will return the proteome with
#' the most proteins that is on uniprot.
#'
#' @param organism_ids NCBI taxa ids
#' @param parallel whether to use parallel processing TRUE or FALSE.
#'
#' @returns a tibble containing the columns Proteome ID, Organism , Organism ID,
#' Protein Count, and reference TRUE or FALSE.
#'
#' @export
#'
#' @examples
get_proteome_ids_from_organism_ids <- function(organism_ids,
                                               parallel = TRUE) {

  # Let the user know if there are duplicate ids
  if(length(organism_ids) - length(unique(organism_ids)) > 0){
    message("There are ", length(organism_ids) - length(unique(organism_ids)),
            " duplicate organism ids in the input. \n")
    message("Only unique ids will be processed. \n")
  }

  organism_ids = unique(organism_ids)

  if (parallel) {
    future::plan(future::multisession, workers = future::availableCores() - 1)
    results <- furrr::future_map_dfr(organism_ids,
                                     get_proteome_id_from_organism_id,
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
  message("Getting Reference Proteomes \n")
  reference <- get_all_reference_proteomes() |> dplyr::pull(Proteome_ID)
  message("Reference Proteomes Sucessfully Retrieved. \n")

  # Annotate with whether it is a reference proteome or not.
  results <- results |>
    dplyr::mutate(reference = dplyr::case_when(
      `Proteome Id` %in% reference ~ TRUE,
      TRUE ~ FALSE
    )) |>
    dplyr::group_by(`Organism Id`) |>
    dplyr::filter(if (any(reference)) reference else `Protein count` == max(`Protein count`)) |>
    dplyr::ungroup()
  # Dealing with missing IDs
  missing_ids = organism_ids[organism_ids %!in% results$`Organism Id`]
  # Print out number that were missing
  if(length(missing_ids > 0)){
  message(length(missing_ids)," organism ids do not have a Uniprot proteome. \n")
  # Print out the ids that they include.
  message("These include the following NCBI ids: ",paste(missing_ids, collapse = ", "))
  missing_ids <- tibble::tibble(`Organism Id` = as.numeric(missing_ids))
  results = dplyr::bind_rows(results,missing_ids)
  }
  return(results)
}
