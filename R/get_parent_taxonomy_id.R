#' Get Parent Taxonomy ID from UniProt Taxonomy API
#'
#' Queries the UniProt taxonomy REST API for a given organism (NCBI taxonomy ID)
#' and returns the parent taxon ID and rank. On failure returns a one-row
#' tibble with NAs and a warning.
#'
#' @param organism_id Integer or character. NCBI taxonomy ID (e.g. 9606 for
#'   human).
#'
#' @return A tibble with columns `parent_id`, `child_id`, `child_rank`. On
#'   success, `parent_id` is the parent taxon ID and `child_rank` is the rank
#'   of the requested organism.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' get_parent_taxonomy_id(9606)   # human
#' get_parent_taxonomy_id(562)    # E. coli
#' }
get_parent_taxonomy_id <- function(organism_id) {
  taxon_url <- paste0("https://rest.uniprot.org/taxonomy/", organism_id)

  resp <- tryCatch(
    httr2::request(taxon_url) |>
      httr2::req_headers(accept = "application/json") |>
      httr2::req_perform(),
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(resp) || httr2::resp_status(resp) != 200) {
    warning(
      "Request failed for ", organism_id,
      if (!is.null(resp)) paste0(" (status: ", httr2::resp_status(resp), ")")
    )
    return(tibble::tibble(
      parent_id = NA_integer_,
      child_id = organism_id,
      child_rank = NA_character_
    ))
  }

  out <- httr2::resp_body_json(resp)

  return(tibble::tibble(
    parent_id = out$parent$taxonId,
    child_id = organism_id,
    child_rank = out$rank
  ))
}
