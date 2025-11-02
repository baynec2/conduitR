#' Title
#'
#' @param proteome_ids
#'
#' @returns
#' @export
#'
#' @examples
get_proteome_taxids_and_types <- function(proteome_ids, parallel = FALSE) {

  # Construct requests
  reqs <- proteome_ids |>
    lapply(function(id) {
      httr2::request(paste0("https://rest.uniprot.org/proteomes/", id))
    })

  # Perform requests (parallel or sequential)
  if (parallel) {
    resps <- httr2::req_perform_parallel(reqs)
  } else {
    resps <- lapply(reqs, httr2::req_perform)
  }

  # Parse each response
  parse_resp <- function(resp, id) {
    if (httr2::resp_status(resp) != 200) {
      warning(paste("Failed to retrieve data for", id))
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

