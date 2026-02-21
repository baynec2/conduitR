#' Get Organism IDs and Proteome Types for UniProt Proteome IDs
#'
#' For each UniProt proteome ID, fetches the NCBI taxonomy ID (organism_id) and
#' proteome type from the UniProt REST API. Used to annotate downloaded
#' proteomes (e.g. in `download_fasta_from_proteome_ids`).
#'
#' @param proteome_ids Character vector of UniProt proteome IDs (e.g.
#'   `"UP000005640"`).
#' @param parallel Logical. If `TRUE`, perform HTTP requests in parallel via
#'   `mapply`. Default `FALSE`.
#'
#' @return A tibble with columns `proteome_id`, `organism_id`, `proteome_type`.
#'   Failed requests have `NA` for `organism_id` and `proteome_type`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' get_proteome_taxids_and_types(c("UP000005640", "UP000000625"))
#' get_proteome_taxids_and_types(c("UP000005640"), parallel = TRUE)
#' }
get_proteome_taxids_and_types <- function(proteome_ids) {

  # Construct requests
  reqs <- proteome_ids |>
    lapply(function(id) {
      httr2::request(paste0("https://rest.uniprot.org/proteomes/", id))
    })

  # Perform requests (parallel or sequential) with tryCatch
  perform_req <- function(req, id) {
    tryCatch({
      resp <- httr2::req_perform(req)
      resp
    }, error = function(e) {
      warning(paste("Request failed for", id, ":", e$message))
      return(NULL)
    })
  }
   resps <- mapply(perform_req, reqs, proteome_ids, SIMPLIFY = FALSE)
  

  # Parse each response
  parse_resp <- function(resp, id) {
    if (is.null(resp) || httr2::resp_status(resp) != 200) {
      return(list(organism_id = NA, proteome_type = NA))
    }
    data <- httr2::resp_body_json(resp)
    taxid <- if (!is.null(data$taxonomy$taxonId)) as.integer(data$taxonomy$taxonId) else NA
    prot_type <- if (!is.null(data$proteomeType)) data$proteomeType else NA
    list(organism_id = taxid, proteome_type = prot_type)
  }

  info_list <- Map(parse_resp, resps, proteome_ids)

  # Combine into tibble
  out <- tibble::tibble(
    proteome_id = proteome_ids,
    organism_id = sapply(info_list, `[[`, "organism_id"),
    proteome_type = sapply(info_list, `[[`, "proteome_type")
  )

  return(out)
}
