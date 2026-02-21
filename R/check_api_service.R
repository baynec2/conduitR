#' Check That an API Endpoint Is Reachable
#'
#' Performs a simple HTTP request to the given URL. If the request succeeds
#' (status 200), the function returns invisibly; otherwise it throws an error.
#' Use this before running workflows that depend on external APIs (e.g. UniProt,
#' NCBI) to fail fast with a clear message when the service is down.
#'
#' @param api_url Character string. URL to test (default: UniProt human proteome
#'   endpoint).
#'
#' @return Invisibly returns `TRUE` on success. Stops with an error if the
#'   request fails or the URL is unreachable.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Default: check UniProt proteomes API
#' check_api_service()
#'
#' # Check a specific UniProt proteome
#' check_api_service("https://rest.uniprot.org/proteomes/UP000005640")
#'
#' # Check NCBI taxonomy API
#' check_api_service(
#'   "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=taxonomy"
#' )
#' }
check_api_service <- function(api_url = "https://rest.uniprot.org/proteomes/UP000005640") {

  # Check API
  api_ok <- tryCatch({
    resp <- httr2::request(api_url) |> httr2::req_perform()
    httr2::resp_status(resp) == 200
  }, error = function(e) {
    FALSE
  })

  if (!api_ok) {
    stop(paste0("API at ", api_url, " is unreachable.
                \nCannot proceed with API-dependent workflows."))
  }

  message(paste0("API at ", api_url, " is avalible. It is okay to proceed."))
  invisible(TRUE)

}
