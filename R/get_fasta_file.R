#' get_fasta_file
#' download proteome file from uniprot.
#' @param proteome_id proteome id
#' @param fasta_fp file path to write to
#'
#' @returns
#' @export
#'
#' @examples
get_fasta_file <- function(proteome_id,
                           fasta_dir = getwd()) {
  # Define the base url
  base_url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?query="
  )

  req <- httr2::request(base_url) |>
    httr2::req_url_query(
      query = paste0("proteome:", proteome_id),
      format = "fasta"
    ) |>
    httr2::req_perform()

  # Define the fasta file path
  fasta_fp <- paste0(fasta_dir, "/", proteome_id, ".fasta")

  # Defining a data frame with information about download

  out <- tibble::tibble(
    proteome_id = proteome_id,
    resp_status = httr2::resp_status(req)
  )

  # Save FASTA content to file
  if (httr2::resp_status(req) < 400 & httr2::resp_has_body(req)) {
    writeLines(httr2::resp_body_string(req), fasta_fp)
    message(
      "Fasta File for Proteome ID: ",
      proteome_id, " sucessfully downloaded"
    )
  } else {
    message("Failed to download FASTA for Proteome: ", proteome_id)
    message("status code: ", httr2::resp_status(req))
    message("has body: ", httr2::resp_has_body(req))
  }

  return(out)
}
