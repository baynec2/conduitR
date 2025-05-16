#' Get Proteome IDs for Multiple Organisms
#'
#' Retrieves proteome IDs from UniProt for a list of NCBI taxonomy IDs, with
#' support for parallel processing and reference proteome prioritization. This
#' function is particularly useful for preparing data for metaproteomics analysis
#' or creating custom protein sequence databases.
#'
#' @param organism_ids A numeric vector of NCBI taxonomy IDs (e.g., 9606 for human,
#'   562 for E. coli). These IDs are used to identify organisms in UniProt.
#' @param parallel Logical indicating whether to use parallel processing (default: TRUE).
#'   When TRUE, uses all available CPU cores minus one for faster processing.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Proteome Id: UniProt proteome identifier
#'     \item Organism: Scientific name of the organism
#'     \item Organism Id: NCBI taxonomy ID
#'     \item Protein count: Number of proteins in the proteome
#'     \item reference: Logical indicating whether this is a reference proteome
#'   }
#'   The function prioritizes reference proteomes when available, and for organisms
#'   without reference proteomes, it selects the proteome with the highest protein
#'   count. Missing organism IDs are included in the output with NA values for
#'   proteome information.
#'
#' @export
#'
#' @examples
#' # Get proteome IDs for human and E. coli:
#' # proteomes <- get_proteome_ids_from_organism_ids(c(9606, 562))
#' 
#' # Process a larger set of organisms in parallel:
#' # all_proteomes <- get_proteome_ids_from_organism_ids(
#' #   c(9606, 562, 10090, 10116),  # Human, E. coli, Mouse, Rat
#' #   parallel = TRUE
#' # )
#' 
#' # Use the results to download FASTA files:
#' # download_fasta_from_organism_ids(proteomes$`Proteome Id`)
#' 
#' # Filter for reference proteomes only:
#' # reference_proteomes <- proteomes[proteomes$reference, ]
#'
#' @note
#' This function:
#' \itemize{
#'   \item Automatically removes duplicate organism IDs
#'   \item Prioritizes reference proteomes when available
#'   \item Falls back to the largest proteome for non-reference organisms
#'   \item Reports missing organism IDs
#'   \item Supports parallel processing for faster execution
#' }
#' 
#' Important considerations:
#' \itemize{
#'   \item The function requires an internet connection to access UniProt
#'   \item API rate limits may affect processing speed
#'   \item Some organisms may not have proteomes in UniProt
#'   \item Reference proteomes are updated every 8 weeks
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[get_proteome_id_from_organism_id]{get_proteome_id_from_organism_id}}
#'     for processing a single organism
#'   \item \code{\link[get_all_reference_proteomes]{get_all_reference_proteomes}}
#'     for retrieving reference proteome information
#'   \item \code{\link[download_fasta_from_organism_ids]{download_fasta_from_organism_ids}}
#'     for downloading proteome FASTA files
#'   \item \code{\link[get_fasta_file]{get_fasta_file}} for downloading
#'     individual FASTA files
#' }
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
