#' Get Proteome ID for a Single Organism
#'
#' Retrieves the proteome ID for a single organism from UniProt using its NCBI
#' taxonomy ID. This function is typically used internally by
#' \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#' for batch processing of multiple organisms.
#'
#' @param id A single NCBI taxonomy ID (e.g., 9606 for human, 562 for E. coli).
#'   This ID is used to query UniProt's proteome database.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item Proteome Id: UniProt proteome identifier
#'     \item Organism: Scientific name of the organism
#'     \item Organism Id: NCBI taxonomy ID
#'     \item Protein count: Number of proteins in the proteome
#'   }
#'   If the request fails, returns NULL and prints an error message.
#'
#' @export
#'
#' @examples
#' # Get proteome ID for human:
#' # human_proteome <- get_proteome_id_from_organism_id(9606)
#' 
#' # Get proteome ID for E. coli:
#' # ecoli_proteome <- get_proteome_id_from_organism_id(562)
#' 
#' # The function is typically used with get_proteome_ids_from_organism_ids:
#' # proteomes <- get_proteome_ids_from_organism_ids(c(9606, 562))
#'
#' @note
#' This function:
#' \itemize{
#'   \item Makes a single API request to UniProt
#'   \item Returns only the first proteome if multiple are found
#'   \item Handles API errors gracefully
#'   \item Is designed for use in batch processing
#' }
#' 
#' For processing multiple organisms, use
#' \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#' which includes additional features like:
#' \itemize{
#'   \item Parallel processing
#'   \item Reference proteome prioritization
#'   \item Duplicate ID handling
#'   \item Missing ID reporting
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[get_proteome_ids_from_organism_ids]{get_proteome_ids_from_organism_ids}}
#'     for processing multiple organisms
#'   \item \code{\link[get_all_reference_proteomes]{get_all_reference_proteomes}}
#'     for retrieving reference proteome information
#'   \item \code{\link[get_fasta_file]{get_fasta_file}} for downloading
#'     proteome FASTA files
#' }
get_proteome_id_from_organism_id <- function(id) {
  # Base URL for request
  base_url <- "https://rest.uniprot.org/proteomes/search?query="
  query <- paste0(paste0("organism_id:", id))
  req <- httr2::request(base_url) |>
    httr2::req_url_query(
      query = query,
      format = "tsv" # Max limit per request
    ) |>
    httr2::req_perform()

  # Parse response
  if (httr2::resp_status(req) == 200) {
    out <- httr2::resp_body_string(req) |>
      readr::read_tsv(show_col_types = FALSE,col_types = "ccnn") |>
      # If there are two proteome ids with the same length, it will return both.
      # Only grabbing first result
      dplyr::slice(1)
    return(out)
  } else {
    message("Request failed for: ", id, "status code: ", httr2::resp_status(req),
            "has body: ", httr2::resp_has_body(req))
    return(NULL)
  }

}
