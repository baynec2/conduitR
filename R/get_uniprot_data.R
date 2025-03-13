#' get_uniprot_data
#'
#' @param ids
#'
#' @returns
#' @export
#'
#' @examples
get_uniprot_data <- function(ids, columns = NULL, batch_size = 150) {
  # Base URL for request
  base_url <- "https://rest.uniprot.org/uniprotkb/search?"
  query <- paste0(paste0("accession:", ids, collapse = " OR "))
  req <- httr2::request(base_url) |>
    httr2::req_url_query(
      query = query,
      format = "tsv",
      fields = columns,
      size = batch_size # Max limit per request
    ) |>
    httr2::req_perform()

  # Parse response
  if (httr2::resp_status(req) == 200) {
    httr2::resp_body_string(req) |>
      readr::read_tsv(show_col_types = FALSE)
  } else {
    message("Request failed, status code: ", httr2::resp_status(req))
    return(NULL)
  }
}
