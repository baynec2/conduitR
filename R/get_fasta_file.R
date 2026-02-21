#' Download Proteome FASTA File from UniProt
#'
#' Downloads a proteome FASTA file from UniProt using their REST API. This function
#' retrieves protein sequences for a specified proteome ID and saves them to a local
#' file. The downloaded file contains protein sequences with UniProt-style headers
#' including protein IDs, organism information, and other metadata.
#'
#' @param proteome_id Character string specifying the UniProt proteome identifier
#'   (e.g., "UP000005640" for human proteome). These IDs can be found in the
#'   UniProt proteomes database.
#' @param fasta_dir Character string specifying the directory where the FASTA file
#'   should be saved. Defaults to the current working directory.
#'
#' @return A tibble containing download information with the following columns:
#'   \itemize{
#'     \item proteome_id: The input proteome identifier
#'     \item resp_status: HTTP response status code (200 for success)
#'     \item source: Either "uniprotkb", "uniparc", or "not_downloaded"
#'   }
#'   The function also saves a FASTA file named "{proteome_id}.fasta" in the
#'   specified directory when sequences are found in UniProtKB or (if empty)
#'   in UniParc.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download human proteome to current directory
#' result <- get_fasta_file("UP000005640")
#'
#' # Download to a specific directory
#' result <- get_fasta_file("UP000005640", fasta_dir = "data/proteomes")
#'
#' # Check download status
#' if (result$resp_status == 200) {
#'   message("Download successful")
#' } else {
#'   message("Download failed")
#' }
#' }
#'
#' @note
#' This function:
#' \itemize{
#'   \item Requires an internet connection to access UniProt
#'   \item Uses the UniProt REST API
#'   \item Creates a new file or overwrites existing files
#'   \item Handles API errors with informative messages
#' }
#' The UniProt API is free to use but has rate limits. For bulk downloads,
#' consider implementing appropriate delays between requests.
get_fasta_file <- function(proteome_id,
                           fasta_dir = getwd()) {
  # Log the proteome id and the time to console.
  log_with_timestamp(paste0("Processing proteome ID: ", proteome_id))

  fasta_fp <- paste0(fasta_dir, "/", proteome_id, ".fasta")

  # 1) Try UniProtKB stream first
  base_url <- "https://rest.uniprot.org/uniprotkb/stream"
  req <- httr2::request(base_url) |>
    httr2::req_url_query(
      query = paste0("proteome:", proteome_id),
      format = "fasta"
    ) |>
    httr2::req_perform()

  status <- httr2::resp_status(req)
  body <- if (httr2::resp_has_body(req)) httr2::resp_body_string(req) else ""

  # Consider "has sequences" if body contains at least one FASTA header
  has_sequences <- nzchar(trimws(body)) && (grepl("^>", body) || grepl("\n>", body, fixed = TRUE))

  if (status < 400L && has_sequences) {
    writeLines(body, fasta_fp)
    log_with_timestamp(paste0(
      "Fasta File for Proteome ID: ",
      proteome_id, " successfully downloaded (UniProtKB)"
    ))
    return(tibble::tibble(
      proteome_id = proteome_id,
      resp_status = status,
      source = "uniprotkb"
    ))
  }

  # 2) If UniProtKB empty or failed, try UniParc proteome stream (better formatted)
  uniparc_url <- paste0(
    "https://rest.uniprot.org/uniparc/proteome/", proteome_id, "/stream"
  )
  req_uniparc <- tryCatch(
    httr2::request(uniparc_url) |>
      httr2::req_url_query(
        compressed = "false",
        format = "fasta"
      ) |>
      httr2::req_perform(),
    error = function(e) NULL
  )

  if (!is.null(req_uniparc) &&
      httr2::resp_status(req_uniparc) < 400L &&
      httr2::resp_has_body(req_uniparc)) {
    body_uniparc <- httr2::resp_body_string(req_uniparc)
    if (nzchar(trimws(body_uniparc)) && (grepl("^>", body_uniparc) || grepl("\n>", body_uniparc, fixed = TRUE))) {
      writeLines(body_uniparc, fasta_fp)
      log_with_timestamp(paste0(
        "Fasta File for Proteome ID: ",
        proteome_id, " successfully downloaded (UniParc)"
      ))
      return(tibble::tibble(
        proteome_id = proteome_id,
        resp_status = httr2::resp_status(req_uniparc),
        source = "uniparc"
      ))
    }
  }

  # 3) Neither source had sequences
  log_with_timestamp(paste0(
    "Failed to download FASTA for Proteome: ", proteome_id
  ))
  message("status code: ", status)
  message("has body: ", httr2::resp_has_body(req))
  tibble::tibble(
    proteome_id = proteome_id,
    resp_status = status,
    source = "not_downloaded"
  )
}
