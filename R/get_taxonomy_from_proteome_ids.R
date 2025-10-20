#' Retrieve NCBI Taxonomy IDs for UniProt Proteome IDs
#'
#' This function takes one or more UniProt proteome IDs and retrieves the
#' corresponding NCBI taxonomy IDs using the UniProt REST API.
#'
#' @param proteome_ids A character vector of UniProt proteome IDs (e.g., "UP000000625").
#' @return A tibble with one row per proteome ID and columns:
#'   \describe{
#'     \item{proteome_id}{The input UniProt proteome ID.}
#'     \item{taxon_id}{The corresponding NCBI taxonomy ID (integer), or NA if not found.}
#'   }
#' @details
#' The function queries the UniProt REST API endpoint
#' `https://rest.uniprot.org/proteomes/{proteome_id}` for each input proteome ID.
#' Invalid IDs or failed requests will return `NA` for `taxon_id`.
#'
#' @examples
#' get_taxonomy_from_proteome_ids("UP000000625")
#' get_taxonomy_from_proteome_ids(c("UP000000625", "UP000005640"))
#'
#' @export
get_taxonomy_from_proteome_ids <- function(proteome_ids) {
  if (missing(proteome_ids) || length(proteome_ids) == 0) {
    stop("Please provide one or more UniProt proteome IDs.")
  }

  proteome_ids <- as.character(proteome_ids)  # ensure character vector

  fetch_one <- function(proteome_id) {
    url <- paste0("https://rest.uniprot.org/proteomes/", proteome_id)

    # safe request
    resp <- tryCatch(
      httr2::request(url) |>
        httr2::req_perform() |>
        httr2::resp_check_status(),
      error = function(e) return(NULL)
    )

    if (is.null(resp)) {
      tibble::tibble(
        proteome_id = proteome_id,
        organism_id = NA_integer_
      )
    } else {
      data <- httr2::resp_body_json(resp)
      tibble::tibble(
        proteome_id = proteome_id,
        organism_id = if (!is.null(data$taxonomy$taxonId)) as.integer(data$taxonomy$taxonId) else NA_integer_
      )
    }
  }

  # loop over all proteome IDs
  result <- purrr::map_dfr(proteome_ids, fetch_one)
  return(result)
}
