#' get_batch_of_proteome_ids
#'
#' get a batch of proteome ids from organism ids
#'
#' Note this can only handle
#'
#' @param ids
#' @param base_url
#'
#' @returns
#' @export
#'
#' @examples
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
