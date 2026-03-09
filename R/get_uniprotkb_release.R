#' Get the current UniProtKB release version
#'
#' Queries the UniProt REST API and returns the current release identifier
#' from the response header. Requires internet access.
#'
#' @return A character string of the form \code{"YYYY_MM"} (e.g. \code{"2026_01"}),
#'   or \code{NULL} if the request fails.
#' @export
get_uniprotkb_release <- function() {
  resp <- httr2::request("https://rest.uniprot.org/uniprotkb/search") |>
    httr2::req_url_query(query = "insulin", size = 1, format = "json") |>
    httr2::req_error(is_error = \(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200) {
    warning("UniProt request failed with status ", httr2::resp_status(resp))
    return(NULL)
  }

  httr2::resp_header(resp, "x-uniprot-release")
}
