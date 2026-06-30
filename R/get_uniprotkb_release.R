#' Get the current UniProtKB release version
#'
#' Queries the UniProt REST API and returns the current release identifier
#' from the response header. Requires internet access.
#'
#' @return A character string of the form \code{"YYYY_MM"} (e.g. \code{"2026_01"}),
#'   or \code{NULL} if the request fails.
#' @export
get_uniprotkb_release <- function() {
  # We only need the `x-uniprot-release` response header, not the search
  # results, so query a single stable accession (P01308, human insulin) rather
  # than a free-text term — deterministic and independent of the search index.
  #
  # Resilient against transient blips: retry on connection-level failures
  # (timeouts, resets) as well as transient HTTP statuses, bound the wait with
  # a timeout, and never throw — this is provenance metadata and must not be
  # able to abort a build. Connection errors return NULL (caller records NA).
  resp <- tryCatch(
    httr2::request("https://rest.uniprot.org/uniprotkb/search") |>
      httr2::req_url_query(query = "accession:P01308", size = 1, format = "json") |>
      httr2::req_error(is_error = \(r) FALSE) |>
      httr2::req_timeout(seconds = 30) |>
      httr2::req_retry(max_tries = 5, retry_on_failure = TRUE) |>
      httr2::req_perform(),
    error = function(e) {
      warning("UniProt release request failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(resp)) {
    return(NULL)
  }

  if (httr2::resp_status(resp) != 200) {
    warning("UniProt request failed with status ", httr2::resp_status(resp))
    return(NULL)
  }

  httr2::resp_header(resp, "x-uniprot-release")
}
