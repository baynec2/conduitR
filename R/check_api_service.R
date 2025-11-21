#' Title
#'
#' @param api_url
#'
#' @returns
#' @export
#'
#' @examples
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
